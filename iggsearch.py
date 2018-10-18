#!/usr/bin/env python

import os, sys, subprocess as sp, numpy as np, time, shutil

class Species:
	def __init__(self):
		self.genes = {}
		self.reads = 0
		self.bases = 0
		self.depth = 0
		self.num_genes = 0
		self.num_covered = 0
		self.fract_covered = 0
		self.length = 0

class Gene:
	def __init__(self):
		self.length = None
		self.reads = 0
		self.bases = 0
		self.depth = 0
		self.weight = None

def parse_arguments():
	""" Parse command line arguments """
	import argparse
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="IGGsearch: metagenome profiling of lineages from the Integrated Human Gut Genomes Database (IGGdb)"
	)
	parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False)
	io = parser.add_argument_group('Input/output options')
	io.add_argument('-o', dest='outdir', type=str, required=True,
		help="""Path to directory to store results.
Directory name should correspond to sample identifier""")
	io.add_argument('-1', type=str, dest='m1', required=True,
		help="""FASTA/FASTQ file containing 1st mate if using paired-end reads.
Otherwise FASTA/FASTQ containing unpaired reads.
Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
Use comma ',' to separate multiple input files""")
	io.add_argument('-2', type=str, dest='m2',
		help="""FASTA/FASTQ file containing 2nd mate if using paired-end reads.""")
	io.add_argument('-d', type=str, dest='db', default=os.environ['IGG_DB'] if 'IGG_DB' in os.environ else None,
		help="""Path to reference database
By default, the MAG_DB environmental variable is used""")

	io.add_argument('-w', dest='weight', choices=['score', 'intra_freq', 'inter_freq', 'none'], default=None,
		help="weight genes based on their score, intra-clade frequency, or inter-clade frequency")
	io.add_argument('--min_intra_freq', type=float, default=0.0)
	io.add_argument('--max_inter_freq', type=float, default=100.0)
	
	io.add_argument('--check', action='store_true', default=False,
		help=argparse.SUPPRESS)
	io.add_argument('--test', action='store_true', default=False,
		help=argparse.SUPPRESS)
	io.add_argument('--max_genes', type=int, help=argparse.SUPPRESS)

	align = parser.add_argument_group('Read alignment options')
	align.add_argument('-n', type=int, dest='max_reads',
		help='# reads to use from input file(s) (use all)')
	align.add_argument('-t', dest='threads', default=1,
		help='Number of threads to use (1)')
	map = parser.add_argument_group('Read filtering options')
	map.add_argument('--mapid', type=float, metavar='FLOAT',
		default=95.0, help='Discard reads with alignment identity < MAPID (95.0)')
	map.add_argument('--aln_cov', type=float, metavar='FLOAT',
		default=0.75, help='Discard reads with alignment coverage < ALN_COV (0.75)')
	args = vars(parser.parse_args())
	return args

def check_args(args):

	# check database
	utility.check_database(args)
	
	if args['test']:
		args['max_reads'] = 1000
		args['max_genes'] = 100000
		args['check'] = True
	
	# create output directory
	if not os.path.isdir(args['outdir']):
		os.makedirs(args['outdir'])
	
	# check input file paths
	for arg in ['m1', 'm2']:
		if not args[arg]:
			continue
		for file in args[arg].split(','):
			if not os.path.isfile(file):
				sys.exit("\nError: Input file does not exist: '%s'\n" % file)

	# input options
	if args['m2'] and not args['m1']:
		sys.exit("\nError: Must specify -1 and -2 if aligning paired end reads\n")

	# sanity check input values
	if args['mapid'] < 1 or args['mapid'] > 100:
		sys.exit("\nError: MAPID must be between 1 and 100\n")
	if args['aln_cov'] < 0 or args['aln_cov'] > 1:
		sys.exit("\nError: ALN_COV must be between 0 and 1\n")

def init_db_info():
	
	db = {}
	file = open('%s/markers.tsv' % args['db'])
	next(file)
	for index, line in enumerate(file):
		
		# truncate database in testing mode
		if args['max_genes'] and index == args['max_genes']: break
		
		# parse line
		row = line.rstrip().split('\t')
		species_id, cluster_id, genome_id, gene_id, rank, length, missing, freq, num_hits, sum_hits, max_hits, all_hits = row
		
		# init objects
		if species_id not in db:
			db[species_id] = Species()

		# TO DO: filter genes
		# ALSO: consider filtering species
		pass

		# store gene
		# TO DO: add more info
		gene = Gene()
		gene.length = int(length)
		gene.intra_freq = float(freq)
		gene.inter_freq = float(sum_hits)
		db[species_id].genes[cluster_id] = gene

	# compute some summary statistics
	for id, sp in db.items():
		sp.num_genes = len(sp.genes)
		if sp.num_genes > 0:
			sp.length = sum([g.length for g in sp.genes.values()])

	print("  total species: %s" % len(db))
	print("  total genes: %s" % sum([len(sp.genes) for sp in db.values()]))
	return db

def weight_genes():

	for species in db.values():

		genes = species.genes.values()

		if len(genes) == 0:
			continue
		
		intra = [g.intra_freq for g in genes]
		intra = [_/sum(intra) for _ in intra]
		
		inter = [g.inter_freq for g in genes]
		if sum(inter) > 0:
			inter = [max(inter) - i for i in inter]
			inter = [i/sum(inter) for i in inter]
		else:
			inter = [1.0/len(genes) for g in genes]
		
		weights = []
		for gene, intra_i, inter_i in zip(genes, intra, inter):
			gene.weight = (intra_i+inter_i)/2.0


def map_reads(args):
	import subprocess
	
	# check output
	out = '%s/mapped_reads.m8' % args['outdir']
	if os.path.exists(out):
		if args['check'] :
			print("  nothing to do")
			return
		else:
			os.remove(out)
	
	# build command
	command = 'python %s/stream_seqs.py' % os.path.dirname(os.path.abspath(sys.argv[0]))
	command += ' -1 %s' % args['m1'] # fasta/fastq
	if args['m2']: command += ' -2 %s' % args['m2'] # mate
	if args['max_reads']: command += ' -n %s' % args['max_reads'] # number of reads
	command += ' | hs-blastn align'
	command += ' -query /dev/stdin'
	command += ' -max_target_seqs 1'
	command += ' -num_alignments 1'
	command += ' -db %s/markers.ffn' % args['db']
	command += ' -outfmt 6'
	command += ' -num_threads %s' % args['threads']
	command += ' -evalue 1e-3'
	command += ' -out %s' % out

	# run command
	if args['verbose']: print("RUNNING: %s" % command)
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	if process.returncode != 0:
		err_message = "\nError encountered executing:\n%s\n\nError message:\n%s\n" % (command, err)
		sys.exit(err_message)
	if args['verbose'] and len(out) > 0: print("hsblastn stdout: %s" % out)
	if args['verbose'] and len(err) > 0: print("hsblastn stderr: %s" % err)

def keep_m8_aln(aln, min_pid, min_aln_cov):
	query_len = aln['query'].split('_')[-1]
	# min pid
	if aln['pid'] < min_pid:
		return False
	# min aln cov
	elif aln['aln']/float(query_len) < min_aln_cov:
		return False
	else:
		return True

def parse_hsblastn(inpath):
	formats = [str,str,float,int,float,float,float,float,float,float,float,float]
	fields = ['query','target','pid','aln','mis','gaps','qstart','qend','tstart','tend','evalue','score']
	for line in open(inpath):
		values = line.rstrip().split()
		yield dict([(field, format(value)) for field, format, value in zip(fields, formats, values)])

def count_reads(args):
	last_query = None
	file = parse_hsblastn('%s/mapped_reads.m8' % args['outdir'])
	for index, aln in enumerate(file):
		
		# skip secondary alignments
		# TO DO: what about paired-end reads?
		if aln['query'] == last_query:
			continue
		last_query = aln['query']
	
		# filter alignments that don't meet cutoffs
		if not keep_m8_aln(aln, args['mapid'], args['aln_cov']):
			continue

		# skip species and/or genes that were excluded from the db
		species_id, gene_id = aln['target'].split('|')
		if species_id not in db or gene_id not in db[species_id].genes:
			continue
		
		# count alignment
		gene = db[species_id].genes[gene_id]
		gene.reads += 1
		gene.bases += aln['aln']
		gene.depth += 1.0*aln['aln']/gene.length

def quantify_abundance():
	total_depth = 0
	for id, sp in db.items():
		genes = sp.genes.values()
		sp.depth = sum([g.depth * g.weight for g in genes])
		sp.reads = sum([g.reads for g in genes])
		sp.fract = sum([g.weight for g in genes if g.reads > 0])
		total_depth += sp.depth
	for id, sp in db.items():
		sp.abun = 100*sp.depth/total_depth if total_depth > 0 else 0.0

def write_profile():
	out = open('%s/species_profile.tsv' % args['outdir'], 'w')
	fields = ['species_id', 'length', 'genes', 'fract', 'reads', 'depth', 'abund']
	out.write('\t'.join(fields)+'\n')
	for spid, sp in db.items():
		row = []
		row.append(spid)
		row.append(sp.length)
		row.append(sp.num_genes)
		row.append(sp.fract)
		row.append(sp.reads)
		row.append(sp.depth)
		row.append(sp.abun)
		out.write('\t'.join([str(_) for _ in row])+'\n')
	out.close()

def log_time(program_start, module_start):
	current_time = time.time()
	program_time = round(current_time - program_start, 2)
	module_time = round(current_time - module_start, 2)
	peak_ram = round(utility.max_mem_usage(), 2)
	print("  module: %s seconds, program: %s seconds, peak RAM: %s GB" % (module_time, program_time, peak_ram))
	
if __name__ == "__main__":

	import time
	program_start = time.time()
	
	args = parse_arguments()

	import utility
	check_args(args) # --> make sure bowtie2 and/or hsblastn is on PATH

	print("\n## Initializing database")
	module_start = time.time()
	db = init_db_info()
	log_time(program_start, module_start)
	
	print("\n## Assigning gene weights")
	module_start = time.time()
	weight_genes()
	log_time(program_start, module_start)

	print("\n## Aligning reads")
	module_start = time.time()
	map_reads(args)
	log_time(program_start, module_start)
	
	print("\n## Counting mapped reads")
	module_start = time.time()
	count_reads(args)
	log_time(program_start, module_start)

	print("\n## Estimating abundance")
	module_start = time.time()
	quantify_abundance()
	log_time(program_start, module_start)

	print("\n## Writing species profile")
	module_start = time.time()
	write_profile()
	log_time(program_start, module_start)




