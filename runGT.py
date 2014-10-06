#!/bin/src/python

import vcf
import updateVcfMongo
import glob
import os

from vcfGT import vcfGT

#directory = r"./files/*.vcf"
directory = r"/data/analysis/avabhyankar/analysis/vcf_for_db/vcf/*.vcf.gz"

for path in glob.glob(directory):
#	vcf_reader = vcf.Reader(open(path,'r'))
	vcf_reader = vcf.Reader(filename=path) 
	ex = vcfGT(vcf_reader)
	ex.InsertIntoGT()

