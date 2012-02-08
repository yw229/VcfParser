#!/usr/bin/env python
'''A VCFv4.0 parser for Python.

The intent of this module is to mimic the ``csv`` module in the Python stdlib,
as opposed to more flexible serialization formats like JSON or YAML.  ``vcf``
will attempt to parse the content of each record based on the data types
specified in the meta-information lines --  specifically the ##INFO and
##FORMAT lines.  If these lines are missing or incomplete, it will check
against the reserved types mentioned in the spec.  Failing that, it will just
return strings.

There main interface is the class: ``Reader``.  It takes a file-like
object and acts as a reader::

    >>> import vcf
    >>> vcf_reader = vcf.Reader(open('test/example.vcf', 'rb'))
    >>> for record in vcf_reader:
    ...     print record
    Record(CHROM=20, POS=14370, REF=G, ALT=['A'])
    Record(CHROM=20, POS=17330, REF=T, ALT=['A'])
    Record(CHROM=20, POS=1110696, REF=A, ALT=['G', 'T'])
    Record(CHROM=20, POS=1230237, REF=T, ALT=[None])
    Record(CHROM=20, POS=1234567, REF=GTCT, ALT=['G', 'GTACT'])


This produces a great deal of information, but it is conveniently accessed.
The attributes of a Record are the 8 fixed fields from the VCF spec::

    * ``Record.CHROM``
    * ``Record.POS``
    * ``Record.ID``
    * ``Record.REF``
    * ``Record.ALT``
    * ``Record.QUAL``
    * ``Record.FILTER``
    * ``Record.INFO``

plus attributes to handle genotype information:

    * ``Record.FORMAT``
    * ``Record.samples``
    * ``Record.genotype``

``samples`` and ``genotype``, not being the title of any column, are left lowercase.  The format
of the fixed fields is from the spec.  Comma-separated lists in the VCF are
converted to lists.  In particular, one-entry VCF lists are converted to
one-entry Python lists (see, e.g., ``Record.ALT``).  Semicolon-delimited lists
of key=value pairs are converted to Python dictionaries, with flags being given
a ``True`` value. Integers and floats are handled exactly as you'd expect::

    >>> vcf_reader = vcf.Reader(open('test/example.vcf', 'rb'))
    >>> record = vcf_reader.next()
    >>> print record.POS
    14370
    >>> print record.ALT
    ['A']
    >>> print record.INFO['AF']
    [0.5]

There are a number of convienience methods and properties for each ``Record`` allowing you to
examine properties of interest::

    >>> print record.num_called, record.call_rate, record.num_unknown
    3 1.0 0
    >>> print record.num_hom_ref, record.num_het, record.num_hom_alt
    1 1 1
    >>> print record.nucl_diversity, record.aaf
    0.6 0.5
    >>> print record.get_hets()
    [Call(sample=NA00002, GT=1|0, GQ=[48])]


``record.FORMAT`` will be a string specifying the format of the genotype
fields.  In case the FORMAT column does not exist, ``record.FORMAT`` is
``None``.  Finally, ``record.samples`` is a list of dictionaries containing the
parsed sample column and ``record.genotype`` is a way of looking up genotypes
by sample name::

    >>> record = vcf_reader.next()
    >>> for sample in record.samples:
    ...     print sample['GT']
    0|0
    0|1
    0/0
    >>> print record.genotype('NA00001')['GT']
    0|0

The genotypes are represented by ``Call`` objects, which have three attributes: the
corresponding Record ``site``, the sample name in ``sample`` and a dictionary of
call data in ``data``::

     >>> call = record.genotype('NA00001')
     >>> print call.site
     Record(CHROM=20, POS=17330, REF=T, ALT=['A'])
     >>> print call.sample
     NA00001
     >>> print call.data
     {'GT': '0|0', 'HQ': [58, 50], 'DP': [3], 'GQ': [49]}

There are also a number of methods::

    >>> print call.called, call.gt_type, call.gt_bases, call.phased
    True 0 T|T True


Metadata regarding the VCF file itself can be investigated through the
following attributes:

    * ``Reader.metadata``
    * ``Reader.infos``
    * ``Reader.filters``
    * ``Reader.formats``
    * ``Reader.samples``

For example::

    >>> vcf_reader.metadata['fileDate']
    '20090805'
    >>> vcf_reader.samples
    ['NA00001', 'NA00002', 'NA00003']
    >>> vcf_reader.filters
    {'q10': Filter(id='q10', desc='Quality below 10'), 's50': Filter(id='s50', desc='Less than 50% of samples have data')}
    >>> vcf_reader.infos['AA'].desc
    'Ancestral Allele'

Random access is supported for files with tabix indexes.  Simply call fetch for the
region you are interested in::

    >>> vcf_reader = vcf.Reader(filename='test/tb.vcf.gz')
    >>> for record in vcf_reader.fetch('20', 1110696, 1230237):
    ...     print record
    Record(CHROM=20, POS=1110696, REF=A, ALT=['G', 'T'])
    Record(CHROM=20, POS=1230237, REF=T, ALT=[None])

Or extract a single row::

    >>> print vcf_reader.fetch('20', 1110696)
    Record(CHROM=20, POS=1110696, REF=A, ALT=['G', 'T'])


The ``Writer`` class provides a way of writing a VCF file.  Currently, you must specify a
template ``Reader`` which provides the metadata::

    >>> vcf_reader = vcf.Reader(filename='test/tb.vcf.gz')
    >>> vcf_writer = vcf.Writer(file('/dev/null', 'w'), vcf_reader)
    >>> for record in vcf_reader:
    ...     vcf_writer.write_record(record)


An extensible script is available to filter vcf files in vcf_filter.py.  VCF filters
declared by other packages will be available for use in this script.  Please
see :doc:`FILTERS` for full description.

'''
import collections
import re
import csv
import gzip
import sys
import itertools

try:
    import pysam
except ImportError:
    pysam = None


VERSION = '0.3.0'


# Metadata parsers/constants
RESERVED_INFO = {
    'AA': 'String', 'AC': 'Integer', 'AF': 'Float', 'AN': 'Integer',
    'BQ': 'Float', 'CIGAR': 'String', 'DB': 'Flag', 'DP': 'Integer',
    'END': 'Integer', 'H2': 'Flag', 'MQ': 'Float', 'MQ0': 'Integer',
    'NS': 'Integer', 'SB': 'String', 'SOMATIC': 'Flag', 'VALIDATED': 'Flag'
}

RESERVED_FORMAT = {
    'GT': 'String', 'DP': 'Integer', 'FT': 'String', 'GL': 'Float',
    'GQ': 'Float', 'HQ': 'Float'
}


_Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc'])
_Filter = collections.namedtuple('Filter', ['id', 'desc'])
_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])


class _vcf_metadata_parser(object):
    '''Parse the metadat in the header of a VCF file.'''
    def __init__(self):
        super(_vcf_metadata_parser, self).__init__()
        self.info_pattern = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>Integer|Float|Flag|Character|String),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.filter_pattern = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.format_pattern = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>.+),
            Description="(?P<desc>.*)"
            >''', re.VERBOSE)
        self.meta_pattern = re.compile(r'''##(?P<key>.+)=(?P<val>.+)''')

    def read_info(self, info_string):
        '''Read a meta-information INFO line.'''
        match = self.info_pattern.match(info_string)
        if not match:
            raise SyntaxError(
                "One of the INFO lines is malformed: %s" % info_string)

        try:
            num = int(match.group('number'))
            if num < 0:
                num = None
        except ValueError:
            num = None

        info = _Info(match.group('id'), num,
                     match.group('type'), match.group('desc'))

        return (match.group('id'), info)

    def read_filter(self, filter_string):
        '''Read a meta-information FILTER line.'''
        match = self.filter_pattern.match(filter_string)
        if not match:
            raise SyntaxError(
                "One of the FILTER lines is malformed: %s" % filter_string)

        filt = _Filter(match.group('id'), match.group('desc'))

        return (match.group('id'), filt)

    def read_format(self, format_string):
        '''Read a meta-information FORMAT line.'''
        match = self.format_pattern.match(format_string)
        if not match:
            raise SyntaxError(
                "One of the FORMAT lines is malformed: %s" % format_string)

        try:
            num = int(match.group('number'))
            if num < 0:
                num = None
        except ValueError:
            num = None

        form = _Format(match.group('id'), num,
                       match.group('type'), match.group('desc'))

        return (match.group('id'), form)

    def read_meta(self, meta_string):
        match = self.meta_pattern.match(meta_string)
        return match.group('key'), match.group('val')


class _Call(object):
    """ A genotype call, a cell entry in a VCF file"""

    def __init__(self, site, sample, data):
        #: The ``_Record`` for this ``_Call``
        self.site = site
        #: The sample name
        self.sample = sample
        #: Dictionary of data from the VCF file
        self.data = data
        self.gt_nums = self.data['GT']
        #: True if the GT is not ./.
        self.called = self.gt_nums is not None

    def __repr__(self):
        return "Call(sample=%s, GT=%s, GQ=%s)" % (self.sample, self.gt_nums, self.data.get('GQ', ''))

    def __eq__(self, other):
        """ Two _Calls are equal if their _Records are equal
            and the samples and ``gt_type``s are the same
        """
        return (self.site == other.site
                and self.sample == other.sample
                and self.gt_type == other.gt_type)

    @property
    def gt_bases(self):
        '''The actual genotype alleles.
           E.g. if VCF genotype is 0/1, return A/G
        '''
        # nothing to do if no genotype call
        if self.called:
            # grab the numeric alleles of the gt string; tokenize by phasing
            phase_char = "/" if not self.phased else "|"
            (a1, a2) = self.gt_nums.split(phase_char)
            # lookup and return the actual DNA alleles
            try:
                return self.site.alleles[int(a1)] + \
                       phase_char + \
                       self.site.alleles[int(a2)]
            except:
                sys.stderr.write("Allele number not found in list of alleles\n")
        else:
            return None

    @property
    def gt_type(self):
        '''The type of genotype.
           hom_ref  = 0
           het      = 1
           hom_alt  = 2  (we don;t track _which+ ALT)
           uncalled = None
        '''
        # extract the numeric alleles of the gt string
        if self.called:
            # grab the numeric alleles of the gt string; tokenize by phasing
            (a1, a2) = self.gt_nums.split("/") \
                if not self.phased else self.gt_nums.split("|")
            # infer genotype type from allele numbers
            if (int(a1) == 0) and (int(a2) == 0): return 0
            elif (int(a1) == 0) and (int(a2) >= 1): return 1
            elif (int(a2) == 0) and (int(a1) >= 1): return 1
            elif (int(a1) >= 1) and (int(a2) >= 1):
                # same alt, so hom_alt
                if a1 == a2: return 2
                # diff alts, so het
                else: return 1
            else: return -1
        else:
            return None

    @property
    def phased(self):
        '''A boolean indicating whether or not
           the genotype is phased for this sample
        '''
        return self.data['GT'] is not None and self.data['GT'].find("|") >= 0

    def __getitem__(self, key):
        """ Lookup value, backwards compatibility """
        return self.data[key]

    @property
    def is_variant(self):
        """ Return True if not a reference call """
        if not self.called:
            return None
        return self.gt_type != 0

    @property
    def is_het(self):
        """ Return True for heterozygous calls """
        if not self.called:
            return None
        return self.gt_type == 1


class _Record(object):
    """ A set of calls at a site.  Equivalent to a row in a VCF file.

        The standard VCF fields CHROM, POS, ID, REF, ALT, QUAL, FILTER,
        INFO and FORMAT are available as properties.

        The list of genotype calls is in the ``samples`` property.
    """

    def __init__(self, CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample_indexes, samples=None):
        self.CHROM = CHROM
        self.POS = POS
        self.ID = ID
        self.REF = REF
        self.ALT = ALT
        self.QUAL = QUAL
        self.FILTER = FILTER
        self.INFO = INFO
        self.FORMAT = FORMAT
        #: 0-based start coordinate
        self.start = self.POS - 1
        #: 1-based end coordinate
        self.end = self.start + len(self.REF)
        #: list of alleles. [0] = REF, [1:] = ALTS
        self.alleles = [self.REF]
        self.alleles.extend(self.ALT)
        #: list of ``_Calls`` for each sample ordered as in source VCF
        self.samples = samples
        self._sample_indexes = sample_indexes

    def __eq__(self, other):
        """ _Records are equal if they describe the same variant (same position, alleles) """
        return (self.CHROM == other.CHROM and
                self.POS == other.POS and
                self.REF == other.REF and
                self.ALT == other.ALT)

    def __iter__(self):
        return iter(self.samples)

    def __str__(self):
        return "Record(CHROM=%(CHROM)s, POS=%(POS)s, REF=%(REF)s, ALT=%(ALT)s)" % self.__dict__

    def add_format(self, fmt):
        self.FORMAT = self.FORMAT + ':' + fmt

    def add_filter(self, flt):
        if self.FILTER is None or self.FILTER == 'PASS':
            self.FILTER = ''
        else:
            self.FILTER = self.FILTER + ';'
        self.FILTER = self.FILTER + flt

    def add_info(self, info, value=True):
        self.INFO[info] = value

    def genotype(self, name):
        """ Lookup a ``_Call`` for the sample given in ``name`` """
        return self.samples[self._sample_indexes[name]]

    @property
    def num_called(self):
        """ The number of called samples"""
        return sum(s.called for s in self.samples)

    @property
    def call_rate(self):
        """ The fraction of genotypes that were actually called. """
        return float(self.num_called) / float(len(self.samples))

    @property
    def num_hom_ref(self):
        """ The number of homozygous for ref allele genotypes"""
        return len([s for s in self.samples if s.gt_type == 0])

    @property
    def num_hom_alt(self):
        """ The number of homozygous for alt allele genotypes"""
        return len([s for s in self.samples if s.gt_type == 2])

    @property
    def num_het(self):
        """ The number of heterozygous genotypes"""
        return len([s for s in self.samples if s.gt_type == 1])

    @property
    def num_unknown(self):
        """ The number of unknown genotypes"""
        return len([s for s in self.samples if s.gt_type is None])

    @property
    def aaf(self):
        """ The allele frequency of the alternate allele.
           NOTE 1: Punt if more than one alternate allele.
           NOTE 2: Denominator calc'ed from _called_ genotypes.
        """
        # skip if more than one alternate allele. assumes bi-allelic
        if len(self.ALT) > 1:
            return None
        hom_ref = self.num_hom_ref
        het = self.num_het
        hom_alt = self.num_hom_alt
        num_chroms = float(2.0*self.num_called)
        return float(het + 2*hom_alt)/float(num_chroms)

    @property
    def nucl_diversity(self):
        """
        pi_hat (estimation of nucleotide diversity) for the site.
        This metric can be summed across multiple sites to compute regional
        nucleotide diversity estimates.  For example, pi_hat for all variants
        in a given gene.

        Derived from:
        \"Population Genetics: A Concise Guide, 2nd ed., p.45\"
          John Gillespie.
        """
        # skip if more than one alternate allele. assumes bi-allelic
        if len(self.ALT) > 1:
            return None
        p = self.aaf
        q = 1.0-p
        num_chroms = float(2.0*self.num_called)
        return float(num_chroms/(num_chroms-1.0)) * (2.0 * p * q)

    def get_hom_refs(self):
        """ The list of hom ref genotypes"""
        return [s for s in self.samples if s.gt_type == 0]

    def get_hom_alts(self):
        """ The list of hom alt genotypes"""
        return [s for s in self.samples if s.gt_type == 2]

    def get_hets(self):
        """ The list of het genotypes"""
        return [s for s in self.samples if s.gt_type == 1]

    def get_unknowns(self):
        """ The list of unknown genotypes"""
        return [s for s in self.samples if s.gt_type is None]


class Reader(object):
    """ Reader for a VCF v 4.0 file, an iterator returning ``_Record objects`` """


    def __init__(self, fsock=None, filename=None, compressed=False, prepend_chr=False):
        """ Create a new Reader for a VCF file.

            You must specify either fsock (stream) or filename.  Gzipped streams
            or files are attempted to be recogized by the file extension, or gzipped
            can be forced with ``compressed=True``
        """
        super(VCFReader, self).__init__()

        if not (fsock or filename):
            raise Exception('You must provide at least fsock or filename')

        if filename:
            self.filename = filename
            if fsock is None:
                self.reader = file(filename)

        if fsock:
            self.reader = fsock
            if filename is None:
                if hasattr(fsock, 'name'):
                    filename = fsock.name
            self.filename = filename

        if compressed or (filename and filename.endswith('.gz')):
            self.reader = gzip.GzipFile(fileobj=self.reader)

        #: metadata fields from header
        self.metadata = None
        #: INFO fields from header
        self.infos = None
        #: FILTER fields from header
        self.filters = None
        #: FORMAT fields from header
        self.formats = None
        self.samples = None
        self._sample_indexes = None
        self._header_lines = []
        self._tabix = None
        self._prepend_chr = prepend_chr
        self._parse_metainfo()

    def __iter__(self):
        return self

    def _parse_metainfo(self):
        '''Parse the information stored in the metainfo of the VCF.

        The end user shouldn't have to use this.  She can access the metainfo
        directly with ``self.metadata``.'''
        for attr in ('metadata', 'infos', 'filters', 'formats'):
            setattr(self, attr, {})

        parser = _vcf_metadata_parser()

        line = self.reader.next()
        while line.startswith('##'):
            self._header_lines.append(line)
            line = line.strip()

            if line.startswith('##INFO'):
                key, val = parser.read_info(line)
                self.infos[key] = val

            elif line.startswith('##FILTER'):
                key, val = parser.read_filter(line)
                self.filters[key] = val

            elif line.startswith('##FORMAT'):
                key, val = parser.read_format(line)
                self.formats[key] = val

            else:
                key, val = parser.read_meta(line.strip())
                self.metadata[key] = val

            line = self.reader.next()

        fields = line.split()
        self.samples = fields[9:]
        self._sample_indexes = dict([(x,i) for (i,x) in enumerate(self.samples)])

    def _map(self, func, iterable, bad='.'):
        '''``map``, but make bad values None.'''
        return [func(x) if x != bad else None
                for x in iterable]

    def _parse_info(self, info_str):
        '''Parse the INFO field of a VCF entry into a dictionary of Python
        types.

        '''
        entries = info_str.split(';')
        retdict = {}
        for entry in entries:
            entry = entry.split('=')
            ID = entry[0]
            try:
                entry_type = self.infos[ID].type
            except KeyError:
                try:
                    entry_type = RESERVED_INFO[ID]
                except KeyError:
                    if entry[1:]:
                        entry_type = 'String'
                    else:
                        entry_type = 'Flag'

            if entry_type == 'Integer':
                vals = entry[1].split(',')
                val = self._map(int, vals)
            elif entry_type == 'Float':
                vals = entry[1].split(',')
                val = self._map(float, vals)
            elif entry_type == 'Flag':
                val = True
            elif entry_type == 'String':
                val = entry[1]

            try:
                if self.infos[ID].num == 1:
                    val = val[0]
            except KeyError:
                pass

            retdict[ID] = val

        return retdict

    def _parse_samples(self, samples, samp_fmt, site):
        '''Parse a sample entry according to the format specified in the FORMAT
        column.'''
        samp_data = []# OrderedDict()
        samp_fmt = samp_fmt.split(':')

        samp_fmt_types = []

        for fmt in samp_fmt:
            try:
                entry_type = self.formats[fmt].type
            except KeyError:
                try:
                    entry_type = RESERVED_FORMAT[fmt]
                except KeyError:
                    entry_type = 'String'
            samp_fmt_types.append(entry_type)

        for name, sample in itertools.izip(self.samples, samples):
            sampdict = self._parse_sample(sample, samp_fmt, samp_fmt_types)
            samp_data.append(_Call(site, name, sampdict))

        return samp_data

    def _parse_sample(self, sample, samp_fmt, samp_fmt_types):
        sampdict = dict([(x, None) for x in samp_fmt])

        for fmt, entry_type, vals in itertools.izip(samp_fmt, samp_fmt_types, sample.split(':')):
            vals = vals.split(',')

            if fmt == 'GT':

                gt = vals[0]

                if gt == './.':
                    sampdict[fmt] = None
                    break
                else:
                    sampdict[fmt] = gt
            else:
                if entry_type == 'Integer':
                    sampdict[fmt] = self._map(int, vals)
                elif entry_type == 'Float' or entry_type == 'Numeric':
                    sampdict[fmt] = self._map(float, vals)
                else:
                    sampdict[fmt] = vals
        return sampdict

    def next(self):
        '''Return the next record in the file.'''
        row = self.reader.next().split()
        chrom = row[0]
        if self._prepend_chr:
            chrom = 'chr' + chrom
        pos = int(row[1])

        if row[2] != '.':
            ID = row[2]
        else:
            ID = None

        ref = row[3]
        alt = self._map(str, row[4].split(','))

        if row[5] == '.':
            qual = None
        else:
            qual = float(row[5]) if '.' in row[5] else int(row[5])
        filt = row[6].split(';') if ';' in row[6] else row[6]
        if filt == 'PASS':
            filt = None
        info = self._parse_info(row[7])

        try:
            fmt = row[8]
        except IndexError:
            fmt = None

        record = _Record(chrom, pos, ID, ref, alt, qual, filt, info, fmt, self._sample_indexes)

        if fmt is not None:
            samples = self._parse_samples(row[9:], fmt, record)
            record.samples = samples

        return record

    def fetch(self, chrom, start, end=None):
        """ fetch records from a Tabix indexed VCF, requires pysam
            if start and end are specified, return iterator over positions
            if end not specified, return individual ``_Call`` at start or None
        """
        if not pysam:
            raise Exception('pysam not available, try "pip install pysam"?')

        if not self.filename:
            raise Exception('Please provide a filename (or a "normal" fsock)')

        if not self._tabix:
            self._tabix = pysam.Tabixfile(self.filename)

        if self._prepend_chr and chrom[:3] == 'chr':
            chrom = chrom[3:]

        # not sure why tabix needs position -1
        start = start - 1

        if end is None:
            self.reader = self._tabix.fetch(chrom, start, start+1)
            try:
                return self.next()
            except StopIteration:
                return None

        self.reader = self._tabix.fetch(chrom, start, end)
        return self


class Writer(object):
    """ VCF Writer """

    fixed_fields = "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT".split()

    def __init__(self, stream, template):
        self.writer = csv.writer(stream, delimiter="\t")
        self.template = template

        for line in template.metadata.items():
            stream.write('##%s=%s\n' % line)
        for line in template.infos.values():
            stream.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % tuple(self._map(str, line)))
        for line in template.formats.values():
            stream.write('##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % tuple(self._map(str, line)))
        for line in template.filters.values():
            stream.write('##FILTER=<ID=%s,Description="%s">\n' % tuple(self._map(str, line)))

        self._write_header()

    def _write_header(self):
        # TODO: write INFO, etc
        self.writer.writerow(self.fixed_fields + self.template.samples)

    def write_record(self, record):
        """ write a record to the file """
        ffs = self._map(str, [record.CHROM, record.POS, record.ID, record.REF]) \
              + [self._format_alt(record.ALT), record.QUAL, record.FILTER or '.',
                 self._format_info(record.INFO), record.FORMAT]

        samples = [self._format_sample(record.FORMAT, sample)
            for sample in record.samples]
        self.writer.writerow(ffs + samples)

    def _format_alt(self, alt):
        return ','.join([x or '.' for x in alt])

    def _format_info(self, info):
        return ';'.join(["%s=%s" % (x, self._stringify(y)) for x, y in info.items()])

    def _format_sample(self, fmt, sample):
        if sample.data["GT"] is None:
            return "./."
        return ':'.join(self._stringify(sample.data[f]) for f in fmt.split(':'))

    def _stringify(self, x, none='.'):
        if type(x) == type([]):
            return ','.join(self._map(str, x, none))
        return str(x) if x is not None else none

    def _map(self, func, iterable, none='.'):
        '''``map``, but make None values none.'''
        return [func(x) if x is not None else none
                for x in iterable]


class Filter(object):
    """ Base class for vcf_filter.py filters """

    name = 'f'
    """ name used to activate filter and in VCF headers """

    description = 'VCF filter base class'
    """ descrtiption used in vcf headers """

    @classmethod
    def customize_parser(self, parser):
        """ hook to extend argparse parser with custom arguments """
        pass

    def __init__(self, args):
        """ create the filter using argparse ``args`` """
        self.threshold = 0

    def __call__(self):
        """ filter a site, return not None if the site should be filtered """
        raise NotImplementedError('Filters must implement this method')


    def filter_name(self):
        """ return the name to put in the VCF header, default is ``name`` + ``threshold`` """
        return '%s%s' % (self.name, self.threshold)


def __update_readme():
    import sys
    file('README.rst', 'w').write(sys.modules[__name__].__doc__)


# backwards compatibility
VCFReader = Reader
VCFWriter = Writer
