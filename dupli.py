#!/bin/src/python
import vcf
import sets 

lasta = []
vcf_reader = vcf.Reader(open('/data/analysis/CasanovaJL/Project_CAS_01015_Exome_2013-09-26/Sample_JL1452/analysis/JL1452.variantCalls.vcf','r'))
from sets import Set
uniq = set()

for r in vcf_reader:
	current =[]
	current.append(r.samples[0].sample)
	current.append(r.CHROM)
	current.append(r.POS)
	current.append(r.REF)
	current = tuple(current)
	if current in uniq:
		print (current,r.ALT,r.REF)
	else:
		uniq.add(current)

#print uniq
