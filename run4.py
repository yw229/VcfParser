import vcf
import updateVcfMongo
import glob
import os

from updateVcfMongo import UpdateVcfMongo

directory = r"./files/*.vcf"
#directory = r"/nethome/avabhyankar/analysis/vcf_for_db/for_yan/*.vcf.gz" 
#directory = r"/data/analysis/avabhyankar/analysis/vcf_for_db/vcf/*.vcf.gz" 

for path in glob.glob(directory):
        vcf_reader = vcf.Reader(open(path,'r'))
	
	#vcf_reader = vcf.Reader(filename=path)

        ex = UpdateVcfMongo(vcf_reader)

        #ex.PrintVcfInsert()

        #ex.PrintVcfUpsert()

        #ex.PrintVcfUpsertFormat()

        ex.PrintVcfUpsertINFO()
