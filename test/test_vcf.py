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


            # check sample name lookup
            for s in self.samples:
                assert site.genotype(s)

            # check ordered access
            self.assertEqual([x.sample for x in site.samples], self.samples)

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

class TestRecord(unittest.TestCase):

    def test_num_calls(self):
        reader = vcf.Reader(fh('example.vcf'))
        for var in reader:
            num_calls = (var.num_hom_ref + var.num_hom_alt + \
                         var.num_het + var.num_unknown)
            self.assertEqual(len(var.samples), num_calls)

    def test_call_rate(self):
        reader = vcf.Reader(fh('example.vcf'))
        for var in reader:
            call_rate = var.call_rate
            if var.POS == 14370:
                self.assertEqual(3.0/3.0, call_rate)
            if var.POS == 17330:
                self.assertEqual(3.0/3.0, call_rate)
            if var.POS == 1110696:
                self.assertEqual(3.0/3.0, call_rate)
            if var.POS == 1230237:
                self.assertEqual(3.0/3.0, call_rate)
            elif var.POS == 1234567:
                self.assertEqual(2.0/3.0, call_rate)

    def test_aaf(self):
        reader = vcf.Reader(fh('example.vcf'))
        for var in reader:
            aaf = var.aaf
            if var.POS == 14370:
                self.assertEqual(3.0/6.0, aaf)
            if var.POS == 17330:
                self.assertEqual(1.0/6.0, aaf)
            if var.POS == 1110696:
                self.assertEqual(None, aaf)
            if var.POS == 1230237:
                self.assertEqual(0.0/6.0, aaf)
            elif var.POS == 1234567:
                self.assertEqual(None, aaf)

    def test_pi(self):
        reader = vcf.Reader(fh('example.vcf'))
        for var in reader:
            pi = var.nucl_diversity
            if var.POS == 14370:
                self.assertEqual(6.0/10.0, pi)
            if var.POS == 17330:
                self.assertEqual(1.0/3.0, pi)
            if var.POS == 1110696:
                self.assertEqual(None, pi)
            if var.POS == 1230237:
                self.assertEqual(0.0/6.0, pi)
            elif var.POS == 1234567:
                self.assertEqual(None, pi)

class TestCall(unittest.TestCase):

    def test_phased(self):
        reader = vcf.Reader(fh('example.vcf'))
        for var in reader:
            phases = [s.phased for s in var.samples]
            if var.POS == 14370:
                self.assertEqual([True, True, False], phases)
            if var.POS == 17330:
                self.assertEqual([True, True, False], phases)
            if var.POS == 1110696:
                self.assertEqual([True, True, False], phases)
            if var.POS == 1230237:
                self.assertEqual([True, True, False], phases)
            elif var.POS == 1234567:
                self.assertEqual([False, False, False], phases)

    def test_gt_bases(self):
        reader = vcf.Reader(fh('example.vcf'))
        for var in reader:
            gt_bases = [s.gt_bases for s in var.samples]
            if var.POS == 14370:
                self.assertEqual(['G|G', 'A|G', 'A/A'], gt_bases)
            elif var.POS == 17330:
                self.assertEqual(['T|T', 'T|A', 'T/T'], gt_bases)
            elif var.POS == 1110696:
                self.assertEqual(['G|T', 'T|G', 'T/T'], gt_bases)
            elif var.POS == 1230237:
                self.assertEqual(['T|T', 'T|T', 'T/T'], gt_bases)
            elif var.POS == 1234567:
                self.assertEqual([None, 'GTCT/GTACT', 'G/G'], gt_bases)

    def test_gt_types(self):
        reader = vcf.Reader(fh('example.vcf'))
        for var in reader:
            gt_types = [s.gt_type for s in var.samples]
            if var.POS == 14370:
                self.assertEqual([0,1,2], gt_types)
            elif var.POS == 17330:
                self.assertEqual([0,1,0], gt_types)
            elif var.POS == 1110696:
                self.assertEqual([1,1,2], gt_types)
            elif var.POS == 1230237:
                self.assertEqual([0,0,0], gt_types)
            elif var.POS == 1234567:
                self.assertEqual([None,1,2], gt_types)

class TestTabix(unittest.TestCase):

    def setUp(self):
        self.reader = vcf.Reader(fh('tb.vcf.gz'))

    def testFetchRange(self):
        lines = list(self.reader.fetch('20', 14370, 14370))
        self.assertEquals(len(lines), 1)
        self.assertEqual(lines[0].POS, 14370)

        lines = list(self.reader.fetch('20', 14370, 17330))
        self.assertEquals(len(lines), 2)
        self.assertEqual(lines[0].POS, 14370)
        self.assertEqual(lines[1].POS, 17330)


        lines = list(self.reader.fetch('20', 1110695, 1234567))
        self.assertEquals(len(lines), 3)

    def testFetchSite(self):
        site = self.reader.fetch('20', 14370)
        assert site.POS == 14370

        site = self.reader.fetch('20', 14369)
        assert site is None




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
        s, out = commands.getstatusoutput('python vcf_filter.py --site-quality 30 test/example.vcf sq')
        #print out
        assert s == 0
        buf = StringIO()
        buf.write(out)
        buf.seek(0)

        print buf.getvalue()
        reader = vcf.Reader(buf)


        # check filter got into output file
        assert 'sq30' in reader.filters

        print reader.filters

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
        s, out = commands.getstatusoutput('python vcf_filter.py --site-quality 30 '
        '--genotype-quality 50 test/example.vcf sq mgq')
        assert s == 0
        #print out
        buf = StringIO()
        buf.write(out)
        buf.seek(0)
        reader = vcf.Reader(buf)

        print reader.filters

        assert 'mgq50' in reader.filters
        assert 'sq30' in reader.filters

class TestRegression(unittest.TestCase):

    def test_issue_16(self):
        reader = vcf.Reader(fh('issue-16.vcf'))
        assert reader.next().QUAL == None



suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestGatkOutput))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFreebayesOutput))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestWriter))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestTabix))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestOpenMethods))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFilter))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(Test1kg))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestRecord))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestCall))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestRegression))
