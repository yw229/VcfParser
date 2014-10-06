#!/bin/src/python
import vcf
import pymongo
import dictMap

class UpdateVcfHeader:
	
	#vcf_reader
	
	def __init__(self,vcf_reader):
		self.vcf_reader = vcf_reader 
            	#global vcf_reader
	
	def ConnectToMongo(self):
		from pymongo import Connection
		connection = Connection()
		db = connection.test
		vcfColl = db.vcfHeader
		return vcfColl
	
	def InsertIntoVcfHeader(self):
		vcf_reader = self.vcf_reader
		vcfColl = self.ConnectToMongo()	
		inf = vcf_reader.infos
		formt = vcf_reader.formats
		filtrs = vcf_reader.filters
		contigs = vcf_reader.contigs
		meta = vcf_reader.metadata
		metaDic = {} 
 
		for i in inf :
			vcfColl.insert({'ID':i,'HeaderTag':'INFO','Type':inf[i].type, 'Number':inf[i].num,'Description':inf[i].desc})

		for i in formt:
			vcfColl.insert({'ID':i,'HeaderTag':'FORMAT','Type':formt[i].type, 'Number':formt[i].num,'Description':formt[i].desc})

		for i in filtrs :
			vcfColl.insert({'ID':i,'HeaderTag':'FILTER','Description':filtrs[i].desc})
		for i in contigs:
			vcfColl.insert({'ID':i,'HeaderTag': 'contig','length':contigs[i].length})
			vcfColl.update({'ID':i,'HeaderTag': 'contig','length':contigs[i].length}, { '$set':{ 'assembly' : 'b37' } },upsert = True, multi =  True )

		for i in meta:
			metaDic[i] = meta[i]
			#vcfColl.update({ 'samples':vcf_reader.samples }, {'$set':metaDic} ,upsert = True, multi =  True )
		
		vcfColl.insert( metaDic )		
		
		return vcfColl


	def PrintVcfInsert(self):
		result = self.InsertIntoVcfHeader()
		for p in result.find():
			print p 	
			
		

	

