#!/usr/bin/python

import sys, gzip

for fname in sys.argv[1:]:
  for line in gzip.open(fname):
      print line,






	
