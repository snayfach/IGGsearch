#!/usr/bin/env python

import os, sys, subprocess as sp, numpy as np, time, shutil, utility

def parse_arguments():
	import argparse
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="IGGsearch: download marker gene database"
	)
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('--out', type=str,
		help="""Directory to download database (current directory)""")
	parser.add_argument('--gut-only', action='store_true', default=False,
		help="""Only download marker genes for 4,558 gut species (False)
Otherwise, download entire database of genes for 23,790 species
The gut-only database is faster and requires less memory""")
	args = vars(parser.parse_args())
	if not args['out']:
		args['out'] = os.getcwd()
	if not os.path.exists(args['out']):
		sys.exit("\nError: specified output directory '%s' does not exist\n" % args['out'])
	return args

def download(args):
	if args['gut_only']:
		url = 'https://www.dropbox.com/s/94ca0jdd3b9naie/iggdb_v1.0.0_gut.tar.gz?dl=0'
		out = '%s/iggdb_v1.0.0_gut.tar.gz' % args['out']
	else:
		url = 'https://www.dropbox.com/s/hewdo8tgq9oyynh/iggdb_v1.0.0.tar.gz?dl=0'
		out = '%s/iggdb_v1.0.0.tar.gz' % args['out']
	import subprocess
	cmd = "wget -O %s %s" % (out, url)
	p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
	out, err = p.communicate()

def main():

	import time, utility
	args = parse_arguments()
	
	program_start = time.time()
	module_start = time.time()
	print("\n## Downloading database")
	download(args)
	utility.log_time(program_start, module_start)





