#!/usr/bin/env python

import os, sys, subprocess as sp, numpy as np, time, shutil

def get_program():
	""" Get program specified by user (species, genes, or snps) """
	if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
		print('Description: Metagenomic species profiling with enhanced coverage of the human gut microbiome')
		print('')
		print('Usage: iggsearch.py <command> [options]')
		print('')
		print('Commands:')
		print('   download download reference database of marker genes')
		print('     search estimate species abundance from a single metagenome')
		print('      merge generate matrix files from multiple runs')
		print('   reformat change the file format of from a single run')
		print('')
		print('Note: use iggsearch.py <command> -h to view usage for a specific command')
		quit()
	elif sys.argv[1] not in ['download', 'search', 'merge', 'reformat']:
		sys.exit("\nError: Unrecognized command: '%s'\n" % sys.argv[1])
		quit()
	else:
		return sys.argv[1]

if __name__ == "__main__":
	
	program = get_program()
	if program == 'search':
		from iggsearch import search
		search.main()
	elif program == 'download':
		from iggsearch import download
		download.main()
	elif program == 'merge':
		sys.exit("\nComing soon...\n")
	elif program == 'reformat':
		sys.exit("\nComing soon...\n")


