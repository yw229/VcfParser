import os
import sys
import vcf

vcf_reader = vcf.VCFReader(open(sys.argv[1]), 'rb')
for var in vcf_reader:
    print var.CHROM, var.POS, var.REF, var.ALT, len(var.samples), \
          var.num_het, var.num_hom_ref, var.num_hom_alt, var.num_unknown
    for s in var.samples:
        if s.called:
            print "\t", s.sample, \
                  s.gt_nums, \
                  s.gt_bases, \
                  s.gt_type, \
                  s.phased
        else:
            print "\t", s.sample, " uncalled"
    print



