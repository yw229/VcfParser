import unittest
import doctest
import os
from StringIO import StringIO

import vcf

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
                assert s in site.genotypes

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

    def testHeader(self):
        print self.reader.samples
        assert self.reader.samples


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


suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestGatkOutput))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFreebayesOutput))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestWriter))
#suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestTabix))

suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestOpenMethods))


