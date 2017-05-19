#!/usr/bin/env python

"""

"""

import os
import sys
import numpy
import argparse


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_dir",
			type=str
		)
	parser.add_argument(
			"out_file",
			type=str
		)
	return parser.parse_args()


def main():
	args = get_args()
	files = list()
	prefiles = os.listdir(args.in_dir)
	for prefile in prefiles:
		if not prefile.startswith('.'):
			files.append(prefile)
	namefile = open("{0}/{1}".format(args.in_dir, files[0]), 'r')
	names = list()
	for line in namefile:
		parts = line.split()
		names.append(parts[0])
	namefile.close()
	outfile = open("{0}".format(args.out_file), 'wb')	
	for name in names:
		rates = list()
		for file in files:
			file = open("{0}/{1}".format(args.in_dir, file), 'r')
			for line in file:
				parts = line.split()
				if parts[0] == name:
					rates.append(float(parts[1]))
			file.close()
		outfile.write("{0}\t{1}\n".format(name, numpy.mean(rates)))
		outfile.flush()
		print name, numpy.mean(rates)
	outfile.close()
					
if __name__ == '__main__':
    main()