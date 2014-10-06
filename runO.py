import glob
import os 
import vcf 
import updateVcfMongo 

from updateVcfMongo import UpdateVcfMongo 

directory = r"./*.vcf"

for path in glob.glob(directory):
	vcf_reader = vcf.Reader(open(path,'r'))
	
	ex = UpdateVcfMongo(vcf_reader)
	
	ex.InsertIntoVcfMongo(5000)
	
	

