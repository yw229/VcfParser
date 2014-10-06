#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import gzip
import sys


def main(infile, outfile):
	if outfile:
		oh = open(outfile, "w")
	else:
		oh = sys.stdout

	fh = gzip.open(infile)
	for line in fh:
		oh.write(line)

	oh.close()


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--file', required=True)
	parser.add_argument('-o', '--outfile')

	args = vars(parser.parse_args())

	main(args['file'], args['outfile'])

