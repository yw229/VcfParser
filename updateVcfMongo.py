#!/bin/src/python
import vcf
import pymongo 
import dictMap

class UpdateVcfMongo:

	def __init__(self,vcf_reader):
		self.vcf_reader = vcf_reader
		#global vcf_reader 
	
	def ConnectToMongo(self):
		from pymongo import Connection
		connection = Connection()
		db = connection.test
		vcfColl = db.VCF
		return vcfColl
		


	def InsertIntoVcfMongo(self,n):
                vcf_reader = self.vcf_reader
               # from dictMap import DictMap
               # dinst = DictMap()
                insertList = []
                vcfColl = self.ConnectToMongo()
                for r in vcf_reader:
                        l = len(r.samples)

                        for i in range(l):
                            insertList.append({'Sample_ID':r.samples[i].sample,'CHROM':r.CHROM,'POS':r.POS,'REF':r.REF,
                                               'alt':str(r.ALT),'QUAL':r.QUAL } )
                            
                            if len(insertList) > n :
					vcfColl.insert(insertList)
                                        insertList = []
                return vcfColl












	
	def UpsertIntoVcfMongo(self):
		vcf_reader = self.vcf_reader
		from dictMap import DictMap
                dinst = DictMap()
                vcfColl = self.ConnectToMongo()
                for r in vcf_reader:
                        l = len(r.samples)
			arr = r.ALT
			f = dinst.my(arr)
			d = dinst.dToString(f)
                        for i in range(l):
				vcfColl.update({'Sample_ID':r.samples[i].sample,'CHROM':r.CHROM,'POS':r.POS,'REF':r.REF,'alt':str(r.ALT) },{'$set':{'QUAL':r.QUAL}}, upsert = True, multi =  True )
				vcfColl.update({'Sample_ID':r.samples[i].sample,'CHROM':r.CHROM,'POS':r.POS,'REF':r.REF,'alt':str(r.ALT) },{'$set': d }, upsert = True, multi =  True )
				                       
                return vcfColl
	
	def UpsertFormatVcfMongo(self):
		vcf_reader = self.vcf_reader
                vcfColl = self.ConnectToMongo()
                for r in vcf_reader:
			FormatKeys = r.FORMAT.split(":")
			sampleSets = r.samples
			l = len(sampleSets)
                        for i in range(l):
				FormatValues = []
				singleSampleSet = sampleSets[i]
				for eachV in singleSampleSet.data:
					FormatValues.append(eachV)
					FormatMap = dict(zip(FormatKeys,FormatValues))	
                       		vcfColl.update({'Sample_ID':singleSampleSet.sample,'CHROM':r.CHROM,'POS':r.POS,'REF':r.REF,'alt':str(r.ALT)},{'$set':FormatMap}, upsert = True, multi =  True )
                return vcfColl

	def UpsertInfoVcfMongo(self):
		vcf_reader = self.vcf_reader
                vcfColl = self.ConnectToMongo()
		for r in vcf_reader:
			infoD = {}
			for eachInfoK in r.INFO:
				#if eachInfoK in ('FS','MQ','MQ0','QD'):
				if eachInfoK in vcf_reader.infos.keys():
					infoD[eachInfoK] = r.INFO.get(eachInfoK)
			sampleSets = r.samples
                        l = len(sampleSets)
			for i in range(l):
				vcfColl.update({'Sample_ID':r.samples[i].sample,'CHROM':r.CHROM,'POS':r.POS,'alt':str(r.ALT),'REF':r.REF },{'$set': infoD }, upsert = True, multi =  True )
		return vcfColl

	def PrintVcfInsert(self):
		result = self.InsertIntoVcfMongo()
#		for p in result.find():
#				print p 
     		
	def PrintVcfUpsert(self):
                result = self.UpsertIntoVcfMongo()
#                for p in result.find():
#                        print p             		
	             		
	def PrintVcfUpsertFormat(self):
		result = self.UpsertFormatVcfMongo()
#		for p in result.find():
#			print p
	
	def PrintVcfUpsertINFO(self):
                result = self.UpsertInfoVcfMongo()
 #               for p in result.find():
 #                       print p             		
