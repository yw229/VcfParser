#!/bin/src/python

import vcf
import pymongo
import dictMap

class GTtest:
	def __init__(self,vcf_reader):
                self.vcf_reader = vcf_reader

	def ConnectToMongo(self):	
		from pymongo import Connection
		conn = Connection()
                db = conn.test
                vcfColl = db.vcfGT
                return vcfColl
        
        def upsertGT(self):
                vcf_reader = self.vcf_reader
                vcfColl = self.ConnectToMongo()
                for r in vcf_reader:
                        for sample in r.samples:
                                vcfColl.update({'CHROM':r.CHROM,'POS':r.POS },
{'$addToSet': { r.genotype(sample.sample)['GT']: sample.sample} },upsert = True, multi =  True)
		return vcfColl 
