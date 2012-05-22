Filtering VCF files
===================

The filter script: vcf_filter.py
--------------------------------

Filtering a VCF file based on some properties of interest is a common enough 
operation that PyVCF offers an extensible script.  ``vcf_filter.py`` does 
the work of reading input, updating the metadata and filtering the records.


Adding a filter
---------------

You can reuse this work by providing a filter class, rather than writing your own filter.
For example, lets say I want to filter each site based on the quality of the site.
I can create a class like this::
   
    import vcf.filters
    class SiteQuality(vcf.filters.Base):

        description = 'Filter sites by quality'
        name = 'sq'

        @classmethod
        def customize_parser(self, parser):
            parser.add_argument('--site-quality', type=int, default=30,
                    help='Filter sites below this quality')

        def __init__(self, args):
            self.threshold = args.site_quality

        def __call__(self, record):
            if record.QUAL < self.threshold:
                return record.QUAL


This class subclasses ``vcf.filters.Base`` which provides the interface for VCF filters.
The ``description``` and ``name`` are metadata about the parser.
The ``customize_parser`` method allows you to add arguments to the script.
We use the ``__init__`` method to grab the argument of interest from the parser.
Finally, the ``__call__`` method processes each record and returns a value if the 
filter failed.  The base class uses the ``name`` and ``threshold`` to create
the filter ID in the VCF file.

To make vcf_filter.py aware of the filter, you can either use the local script option
or declare an entry point.  To use a local script, simply call vcf_filter::

    $ vcf_filter.py --local-script my_filters.py ...

To use an entry point, you need to declare a ``vcf.filters`` entry point in your ``setup``::

    setup(
        ...
        entry_points = {
            'vcf.filters': [
                'site_quality = module.path:SiteQuality',
            ]
        }
    )

Either way, when you call vcf_filter.py, you should see your filter in the list of available filters::

    usage: vcf_filter.py [-h] [--no-short-circuit] [--no-filtered] 
                  [--output OUTPUT] [--local-script LOCAL_SCRIPT]
                  input filter [filter_args] [filter [filter_args]] ...
                

    Filter a VCF file

    positional arguments:
      input                 File to process (use - for STDIN) (default: None)

    optional arguments:
      -h, --help            Show this help message and exit. (default: False)
      --no-short-circuit    Do not stop filter processing on a site if any filter
                            is triggered (default: False)
      --output OUTPUT       Filename to output [STDOUT] (default: <open file
                            '<stdout>', mode 'w' at 0x1002841e0>)
      --no-filtered         Output only sites passing the filters (default: False)
      --local-script LOCAL_SCRIPT
                            Python file in current working directory with the
                            filter classes (default: None)

    sq:
      Filter sites by quality

      --site-quality SITE_QUALITY
                            Filter sites below this quality (default: 30)

The filter base class: vcf.filters.Base
---------------------------------------

.. autoclass:: vcf.filters.Base
   :members:

