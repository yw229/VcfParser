import vcf
import cProfile

def parse_1kg():
    for line in vcf.Reader(filename='1kg.vcf.gz'):
        pass


cProfile.run('parse_1kg()', '1kg.prof')

import pstats
p = pstats.Stats('1kg.prof')

p.strip_dirs().sort_stats('time').print_stats()

