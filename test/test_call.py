import os
import sys
import vcf

vcf_reader = vcf.VCFReader(open(sys.argv[1]), 'rb')
for var in vcf_reader:
    print var.CHROM, var.POS, var.REF, var.ALT
    for sample in vcf_reader.samples:
        if var.genotypes[sample].called:
            print "\t", sample, \
                  var.genotypes[sample].gt_nums, \
                  var.genotypes[sample].gt_bases, \
                  var.genotypes[sample].gt_type, \
                  var.genotypes[sample].phased
        else:
            print "\t", sample, " uncalled"
    print



