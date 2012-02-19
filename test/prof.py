import vcf
import cProfile
import timeit
import pstats
import sys

def parse_1kg():
    for line in vcf.Reader(filename='test/1kg.vcf.gz'):
        pass

if len(sys.argv) == 1:
    sys.argv.append(None)

if sys.argv[1] == 'profile':
    cProfile.run('parse_1kg()', '1kg.prof')
    p = pstats.Stats('1kg.prof')
    p.strip_dirs().sort_stats('time').print_stats()

elif sys.argv[1] == 'time':
    t = timeit.timeit('parse_1kg()',  "from __main__ import parse_1kg", number=1)
    print t
else:
    print 'prof.py profile/time'
