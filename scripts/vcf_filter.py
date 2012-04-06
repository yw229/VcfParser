#!/usr/bin/env python
import sys
import argparse
import pkg_resources

import vcf
from vcf.parser import _Filter

def create_core_parser():
    parser = argparse.ArgumentParser(description='Filter a VCF file',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            add_help=False
            )
    parser.add_argument('-h', '--help', action='store_true',
            help='Show this help message and exit.')
    parser.add_argument('input', metavar='input', type=argparse.FileType('rb'), nargs='?', default=None,
            help='File to process (use - for STDIN)')
    parser.add_argument('filters', metavar='filter', type=str, nargs='*', default=None,
            help='Filters to use')
    parser.add_argument('--no-short-circuit', action='store_true',
            help='Do not stop filter processing on a site if a single filter fails')
    parser.add_argument('--output', action='store', default=sys.stdout,
            help='Filename to output [STDOUT]')
    parser.add_argument('--no-filtered', action='store_true',
            help='Output only sites passing the filters')
    parser.add_argument('--local-script', action='store', default=None,
            help='Python file in current working directory with the filter classes')
    parser.add_argument('rest', nargs=argparse.REMAINDER, help=argparse.SUPPRESS)
    
    return parser

# argument parsing strategy
# loading a script given at the command line poses a difficulty 
# for using the argparse in a simple way -- the command line arguments
# are not completely known the first time command line is parsed
# requirements
#  - display all options grouped by the filters in help screen
#  - check if only arguments for currently used filters are given
#  - to increase legibility when using more filters, arguments should 
#    follow the filter name
#  - it is good to specify the filters explicitly by name, 
#    because the order of filtering can matter
# solution
# - change the command syntax to 
#   vcf_filter.py --core-options input filter1 --filter1-args filter2 filter3
# - parse the core program options with parse_known_args
# - use add_argument_group for filters (subparsers won't work, they require 
#   the second command in argv[1])
# - create all-filters parser when displaying the help
# - parse the arguments incrementally on argparse.REMAINDER of the previous

    # TODO: allow filter specification by short name
    # TODO: flag that writes filter output into INFO column
    # TODO: argument use implies filter use
    # TODO: parallelize
    # TODO: prevent plugins raising an exception from crashing the script

def main():
    # dynamically build the list of available filters
    filters = {}

    # parse command line args
    # (mainly because of local_script)
    parser = create_core_parser()
    (args, unknown_args) = parser.parse_known_args()

    def addfilt(filt):
        filters[filt.name] = filt
        arg_group = parser.add_argument_group(filt.name, filt.description)
        filt.customize_parser(arg_group)
    
    # look for global extensions
    for p in pkg_resources.iter_entry_points('vcf.filters'):
        filt = p.load()
        addfilt(filt)

    # add all classes from local script, if present
    if args.local_script != None:
        import inspect
        import sys, os
        sys.path.insert(0, os.getcwd())
        module_name = args.local_script.replace('.py', '')
        mod = __import__(module_name)
        classes = inspect.getmembers(mod, inspect.isclass)
        for name, cls in classes:
            addfilt(cls)

    # consume all the filter arguments to correctly determine
    # the input argument
    args = parser.parse_args()
    
    # show help after loading of the local module
    # so it can include the additional options
    if args.help or len(args.filters) == 0 or args.input == None:
        parser.print_help()
        parser.exit()

    inp = vcf.Reader(args.input)

    # create an extra parser that keeps only the arguments
    # of selected filters
    filter_parser = argparse.ArgumentParser(description='Agruments of the selected filters', add_help=False)
    for name in args.filters: filters[name].customize_parser(filter_parser)

    # parse the arguments not known to the global arg parser
    # (this raises error if there is some unknown argument)
    print "unk:", unknown_args
    filter_args = filter_parser.parse_args(unknown_args)

    # build filter chain
    chain = []
    for name in args.filters:
        f = filters[name](filter_args)
        chain.append(f)
        inp.filters[f.filter_name()] = _Filter(f.filter_name(), f.description)

    oup = vcf.Writer(args.output, inp)

    # apply filters
    short_circuit = not args.no_short_circuit
    drop_filtered = args.no_filtered

    for record in inp:
        output_record = True
        for filt in chain:
            result = filt(record)
            if result == None: continue

            # save some work by skipping the rest of the code
            if drop_filtered: 
                output_record = False
                break
            
            record.add_filter(filt.filter_name())
            if short_circuit: break

        if output_record:
            # use PASS only if other filter names appear in the FILTER column 
            #FIXME: is this good idea?
            if record.FILTER == '.' and not drop_filtered: record.FILTER = 'PASS'
            oup.write_record(record)

if __name__ == '__main__': main()
