#!/src/bin

import vcf

path = '/nethome/ywang/01-0055.filtered.variantCalls.snpEff.vcf.gz' 
vcf_reader = vcf.Reader(filename=path)
for r in vcf_reader:
	print vcf_reader.next()
