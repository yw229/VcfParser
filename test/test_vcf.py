import unittest
import doctest
import os
import commands
from StringIO import StringIO

import vcf
import vcf_filter


suite = doctest.DocTestSuite(vcf)


def fh(fname):
    return file(os.path.join(os.path.dirname(__file__), fname))

class TestGatkOutput(unittest.TestCase):

    filename = 'gatk.vcf'

    samples = ['BLANK', 'NA12878', 'NA12891', 'NA12892',
            'NA19238', 'NA19239', 'NA19240']
    formats = ['AD', 'DP', 'GQ', 'GT', 'PL']
    infos = ['AC', 'AF', 'AN', 'BaseQRankSum', 'DB', 'DP', 'DS',
            'Dels', 'FS', 'HRun', 'HaplotypeScore', 'InbreedingCoeff',
            'MQ', 'MQ0', 'MQRankSum', 'QD', 'ReadPosRankSum']

    n_calls = 37

    def setUp(self):
        self.reader = vcf.Reader(fh(self.filename))

    def testSamples(self):
        self.assertEqual(self.reader.samples, self.samples)

    def testFormats(self):
        self.assertEqual(set(self.reader.formats), set(self.formats))

    def testInfos(self):
        self.assertEqual(set(self.reader.infos), set(self.infos))


    def testCalls(self):
        n = 0

        for site in self.reader:
            n += 1
            self.assertEqual(len(site.samples), len(self.samples))

            # check sample ordering is preserved:
            self.assertEqual([x['name'] for x in site.samples], self.samples)

            # check sample name lookup
            for s in self.samples:
                assert site.genotype(s)

        self.assertEqual(n,  self.n_calls)


class TestFreebayesOutput(TestGatkOutput):

    filename = 'freebayes.vcf'
    formats = ['AO', 'DP', 'GL', 'GLE', 'GQ', 'GT', 'QA', 'QR', 'RO']
    infos = ['AB', 'ABP', 'AC', 'AF', 'AN', 'AO', 'BVAR', 'CIGAR',
            'DB', 'DP', 'DPRA', 'EPP', 'EPPR', 'HWE', 'LEN', 'MEANALT',
            'NUMALT', 'RPP', 'MQMR', 'ODDS', 'MQM', 'PAIREDR', 'PAIRED',
            'SAP', 'XRM', 'RO', 'REPEAT', 'XRI', 'XAS', 'XAI', 'SRP',
            'XAM', 'XRS', 'RPPR', 'NS', 'RUN', 'CpG', 'TYPE']
    n_calls = 104


class Test1kg(unittest.TestCase):

    def testParse(self):
        reader = vcf.Reader(fh('1kg.vcf.gz'))

        self.assertEqual(len(reader.samples), 629)
        for _ in reader:
            pass


class TestWriter(unittest.TestCase):

    def testWrite(self):

        reader = vcf.Reader(fh('gatk.vcf'))
        out = StringIO()
        writer = vcf.Writer(out, reader)

        records = list(reader)

        map(writer.write_record, records)
        out.seek(0)
        reader2 = vcf.Reader(out)

        self.assertEquals(reader.samples, reader2.samples)
        self.assertEquals(reader.formats, reader2.formats)
        self.assertEquals(reader.infos, reader2.infos)

        for l, r in zip(records, reader2):
            self.assertEquals(l.samples, r.samples)



class TestTabix(unittest.TestCase):

    def setUp(self):
        self.reader = vcf.Reader(fh('tb.vcf.gz'))

    def testFetch(self):
        lines = list(self.reader.fetch('20', 14370-1, 14370+1))
        self.assertEquals(len(lines), 1)
        self.assertEqual(lines[0].POS, 14370)

        lines = list(self.reader.fetch('20', 14370-1, 17330+1))
        self.assertEquals(len(lines), 2)
        self.assertEqual(lines[0].POS, 14370)
        self.assertEqual(lines[1].POS, 17330)


        lines = list(self.reader.fetch('20', 1110695, 1234567))
        self.assertEquals(len(lines), 3)


class TestOpenMethods(unittest.TestCase):

    samples = 'NA00001 NA00002 NA00003'.split()

    def testOpenFilehandle(self):
        r = vcf.Reader(fh('example.vcf'))
        self.assertEqual(self.samples, r.samples)
        self.assertEqual('example.vcf', os.path.split(r.filename)[1])

    def testOpenFilename(self):
        r = vcf.Reader(filename='test/example.vcf')
        self.assertEqual(self.samples, r.samples)

    def testOpenFilehandleGzipped(self):
        r = vcf.Reader(fh('tb.vcf.gz'))
        self.assertEqual(self.samples, r.samples)

    def testOpenFilenameGzipped(self):
        r = vcf.Reader(filename='test/tb.vcf.gz')
        self.assertEqual(self.samples, r.samples)


class TestFilter(unittest.TestCase):


    def testApplyFilter(self):
        out = commands.getoutput('python vcf_filter.py --site-quality 30 test/example.vcf site_quality')
        #print out
        buf = StringIO()
        buf.write(out)
        buf.seek(0)
        reader = vcf.Reader(buf)

        # check filter got into output file
        assert 'sq30' in reader.filters

        # check sites were filtered
        n = 0
        for r in reader:
            if r.QUAL < 30:
                assert 'sq30' in r.FILTER
                n += 1
            else:
                assert 'sq30' not in r.FILTER
        assert n == 2


    def testApplyMultipleFilters(self):
        out = commands.getoutput('python vcf_filter.py --site-quality 30 '
        '--genotype-quality 50 test/example.vcf site_quality min_genotype_quality')
        #print out
        buf = StringIO()
        buf.write(out)
        buf.seek(0)
        reader = vcf.Reader(buf)

        assert 'mgq50' in reader.filters
        assert 'sq30' in reader.filters



suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestGatkOutput))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFreebayesOutput))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestWriter))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestTabix))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestOpenMethods))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFilter))


suite.addTests(unittest.TestLoader().loadTestsFromTestCase(Test1kg))

