#!/usr/bin/env python
import sys
import argparse
import pkg_resources

import vcf
from vcf.parser import _Filter

parser = argparse.ArgumentParser(description='Filter a VCF file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False
        )
parser.add_argument('-h', '--help', action='store_true',
        help='Show this help message and exit.')
parser.add_argument('input', metavar='input', type=str, nargs='?', default=None,
        help='File to process (use - for STDIN)')
parser.add_argument('filters', metavar='filter', type=str, nargs='*', default=None,
        help='Filters to use')
parser.add_argument('--no-short-circuit', action='store_true',
        help='Do not stop filter processing on a site if a single filter fails.')
parser.add_argument('--output', action='store', default=sys.stdout,
        help='Filename to output (default stdout)')
parser.add_argument('--no-filtered', action='store_true',
        help='Remove failed sites')
parser.add_argument('--local-script', action='store', default=None,
        help='File in current working directory with filter classes (with .py)')


    # TODO: allow filter specification by short name
    # TODO: flag that writes filter output into INFO column
    # TODO: argument use implies filter use
    # TODO: parallelize
    # TODO: prevent plugins raising an exception from crashing the script

def main():
    # dynamically build the list of available filters
    filters = {}
    filter_help = '\n\navailable filters:'

    def addfilt(filt, help):
        filters[filt.name] = filt
        filt.customize_parser(parser)
        help += '\n  %s:\t%s' % (filt.name, filt.description)
    
    for p in pkg_resources.iter_entry_points('vcf.filters'):
        filt = p.load()
        addfilt(filt, filter_help)

    # parse command line args
    args = parser.parse_args()

    # add all classes from local script, if present
    if args.local_script != None:
        import inspect
        module_name = args.local_script.replace('.py', '')
        mod = __import__(module_name)
        classes = inspect.getmembers(mod, inspect.isclass)
        for name, cls in classes:
            addfilt(cls, filter_help)

    parser.description += filter_help
    
    # show help after loading of the local module
    # so it can include the additional options
    if args.help or len(args.filters) == 0 or args.input == None:
        parser.print_help()
        parser.exit()
        
    inp = vcf.Reader(file(args.input[0]))

    # build filter chain
    chain = []
    for name in args.filters:
        f = filters[name](args)
        chain.append(f)
        inp.filters[f.filter_name()] = _Filter(f.filter_name(), f.description)

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

        if (not args.no_filtered) or (record.FILTER == '.'):
            oup.write_record(record)

if __name__ == '__main__': main()
