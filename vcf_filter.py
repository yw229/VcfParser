#!/usr/bin/env python
import sys
import argparse
import pkg_resources

import vcf


parser = argparse.ArgumentParser(description='Filter a VCF file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )
parser.add_argument('input', metavar='input', type=str, nargs=1,
        help='File to process (use - for STDIN)')
parser.add_argument('filters', metavar='filter', type=str, nargs='+',
        help='Filters to use')
parser.add_argument('--no-short-circuit', action='store_true',
        help='Do not stop filter processing on a site if a single filter fails.')
parser.add_argument('--output', action='store', default=sys.stdout,
        help='Filename to output (default stdout)')



class SiteQuality(vcf.Filter):

    description = 'Filter sites by quality'
    name = 'site_quality'
    short_name = 'sq'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--site-quality', type=int, default=30,
                help='Filter sites below this quality')

    def __init__(self, args):
        self.threshold = args.site_quality

    def __call__(self, record):
        if record.QUAL < self.threshold:
            return record.QUAL



class VariantGenotypeQuality(vcf.Filter):

    description = 'Demand a minimum quality associated with a non reference call'
    name = 'min_genotype_quality'
    short_name = 'mgq'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--genotype-quality', type=int, default=50,
                help='Filter sites with no genotypes above this quality')

    def __init__(self, args):
        self.threshold = args.genotype_quality

    def __call__(self, record):
        variants = (x for x in record.samples if x['GT'] != '0/0' and x['GT'] != './.')

        vgq = max([max(x['GQ']) for x in variants])
        if vgq < self.threshold:
            return vgq



if __name__ == '__main__':
    # TODO: allow filter specification by short name
    # TODO: flag that writes filter output into INFO column
    # TODO: argument use implies filter use
    # TODO: parallelize
    # TODO: prevent plugins raising an exception from crashing the script


    # dynamically build the list of available filters
    filters = {}
    filter_help = '\n\navailable filters:'

    for p in pkg_resources.iter_entry_points('vcf.filters'):
        filt = p.load()
        filters[filt.name] = filt
        filt.customize_parser(parser)
        filter_help += '\n  %s:\t%s' % (filt.name, filt.description)

    parser.description += filter_help

    # parse command line args
    args = parser.parse_args()

    inp = vcf.Reader(file(args.input[0]))

    # build filter chain
    chain = []
    for name in args.filters:
        f = filters[name](args)
        chain.append(f)
        inp.filters[f.filter_name()] = vcf._Filter(f.filter_name(), f.description)


    oup = vcf.Writer(args.output, inp)

    # apply filters
    short_circuit = not args.no_short_circuit

    for record in inp:
        for filt in chain:
            result = filt(record)
            if result:
                record.add_filter(filt.filter_name())
                if short_circuit:
                    break
        oup.write_record(record)


