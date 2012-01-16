import unittest
import os

import vcf


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
        self.reader = vcf.VCFReader(fh(self.filename))

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
            self.assertEqual([x['name'] for x in site.samples.values()], self.samples)

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



if __name__ == '__main__':
    unittest.main()

