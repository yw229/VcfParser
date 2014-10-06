#!/bin/src/python

import vcf
import glob
import os

from testGT import GTtest

directory = r"./files/*.vcf"

for path in glob.glob(directory):
	vcf_reader = vcf.Reader(open(path,'r'))
	ex = GTtest(vcf_reader)
	ex.upsertGT()

	
