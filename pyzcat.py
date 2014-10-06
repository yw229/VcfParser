#!/usr/bin/python

import sys

def gziplines(fname):
	from subprocess import Popen, PIPE
	f = Popen(['zcat',fname], stdout=PIPE)
	for line in f.stdout:
		yield line 

for fname in sys.argv[1:]:
	for line in gziplines(fname):
		print line,

