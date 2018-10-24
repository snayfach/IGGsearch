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
		help="""Path to reference database. By default, the IGG_DB environmental variable is used""")
		
	stats = parser.add_argument_group('Quantification options')
	stats.add_argument('--presabs_cutoff', type=float, default=0.5,
		help="""Cutoff for determining species presence/absence. By default a species is called present if at least 50% of its marker genes recruit at least 1 mapped read (default=0.5)""")
	stats.add_argument('--min_intra', type=float, default=0.0,
		help="Exclude database genes found in fewer than <min_intra_freq> of reference genomes within species (default=0, range=[0,100])")
	stats.add_argument('--max_inter', type=float, default=100.0,
		help="Exclude database genes found in greater than <max_inter_freq> of reference genomes within species (default=1.0, range=[0,1])")

	speed = parser.add_argument_group('Pipeline speed')
	speed.add_argument('-n', type=int, dest='max_reads',
		help='# reads to use from input file(s) (use all)')
	speed.add_argument('-t', dest='threads', default=1,
		help='Number of threads to use (1)')
	speed.add_argument('--check', action='store_true', default=False, help=argparse.SUPPRESS)
	speed.add_argument('--test', action='store_true', default=False, help=argparse.SUPPRESS)
	speed.add_argument('--max_genes', type=int, help=argparse.SUPPRESS)

	map = parser.add_argument_group('Read filtering options')
	map.add_argument('--mapid', type=float, metavar='FLOAT',
		default=95.0, help='Discard reads with alignment identity < MAPID (95.0)')
	map.add_argument('--aln_cov', type=float, metavar='FLOAT',
		default=0.75, help='Discard reads with alignment coverage < ALN_COV (0.75)')
	map.add_argument('--readq', type=float, metavar='FLOAT',
		default=20.0, help='Minimum average-base-quality per read (20.0)')
	map.add_argument('--mapq', type=float, metavar='FLOAT',
		default=0, help='Minimum map quality score per read (0)')
		
	args = vars(parser.parse_args())
	return args

def which(program):
	""" Mimics unix 'which' function """
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
	for path in os.environ["PATH"].split(os.pathsep):
		path = path.strip('"')
		exe_file = os.path.join(path, program)
		if is_exe(exe_file):
			return True
	return False

def check_database(args):
	if args['db'] is None:
		error = "\nError: No reference database specified\n"
		error += "Use the flag -d to specify a database,\n"
		error += "Or set the IGG_DB environmental variable: export IGG_DB=/path/to/igg_db\n"
		sys.exit(error)
	if not os.path.exists(args['db']):
		error = "\nError: Specified reference database does not exist: %s" % args['db']
		error += "\nCheck that you've entered the path correctly and the database exists\n"
		sys.exit(error)
	for file in []:
		path = '%s/%s' % (args['db'], file)
		if not os.path.exists(path):
			error = "\nError: Could not locate required database file: %s\n" % path
			sys.exit(error)

def check_args(args):

	# check executables
	for exe in ['bowtie2', 'samtools']:
		if not which(exe):
			sys.exit("\nError: required program '%s' not executable or not found on $PATH\n" % exe)

	# check database
	check_database(args)
	
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
	
	num_skipped = 0
	db = {}
	file = open('%s/markers.tsv' % args['db'])
	next(file)
	for index, line in enumerate(file):
		
		# truncate database in testing mode
		if args['max_genes'] and index == args['max_genes']: break
		
		# parse line
		row = line.rstrip().split('\t')
		species_id, cluster_id, length, missing, intra_freq, inter_count, inter_max, inter_sum, inter_freqs = row
		
		# filter genes
		if (float(intra_freq) < args['min_intra'] or float(inter_sum) > args['max_inter']):
			num_skipped += 1
			continue
		
		# TO DO: filter species by genome quality?
		
		# init objects
		if species_id not in db:
			db[species_id] = Species()
		
		# store gene
		# TO DO: add more info
		gene = Gene()
		gene.length = int(length)
		db[species_id].genes[cluster_id] = gene

	# compute some summary statistics
	for id, sp in db.items():
		sp.num_genes = len(sp.genes)
		if sp.num_genes > 0:
			sp.length = sum([g.length for g in sp.genes.values()])

	print("  total species: %s" % len(db))
	print("  total genes: %s" % sum([len(sp.genes) for sp in db.values()]))
	print("  excluded genes: %s" % num_skipped)
	
	return db


def map_reads_bt2(args):
	import subprocess
	
	out = '%s/mapped_reads.bam' % args['outdir']
	
	# Run bowtie2
	command = 'bowtie2 --no-unal '
	command += '-x %s/markers ' % args['db']
	if args['max_reads']: command += '-u %s ' % args['max_reads']
	command += '--threads %s ' % args['threads']
	if args['m2']:
		command += '-1 %s -2 %s ' % (args['m1'], args['m2'])
	else:
		command += '-U %s ' % args['m1']
	# Output sorted bam file
	command += '| samtools view -b - --threads %s ' % args['threads']
	command += '| samtools sort - --threads %s ' % args['threads']
	command += '-O BAM > %s ' % out

	# run command
	print("  running: %s" % command)
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	if process.returncode != 0:
		err_message = "\nError encountered executing:\n%s\n\nError message:\n%s\n" % (command, err)
		sys.exit(err_message)

def keep_aln(aln, min_pid, min_readq, min_mapq, min_aln_cov):
	align_len = len(aln.query_alignment_sequence)
	query_len = aln.query_length
	# min pid
	if 100*(align_len-dict(aln.tags)['NM'])/float(align_len) < min_pid:
		return False
	# min read quality
	elif np.mean(aln.query_qualities) < min_readq:
		return False
	# min map quality
	elif aln.mapping_quality < min_mapq:
		return False
	# min aln cov
	elif align_len/float(query_len)  < min_aln_cov:
		return False
	else:
		return True

def count_reads_bt2(args):
	aligned = 0
	mapped = 0
	import pysam
	bam_path =  '%s/mapped_reads.bam' % args['outdir']
	bamfile = pysam.AlignmentFile(bam_path, "r")
	for index, aln in enumerate(bamfile.fetch(until_eof = True)):
	
		# skip species and/or genes that were excluded from the db
		species_id, gene_id = bamfile.getrname(aln.reference_id).split('|')
		if species_id not in db or gene_id not in db[species_id].genes:
			continue
		
		# filter alignments
		aligned += 1
		args['readq'] = 0
		args['mapq'] = 0
		if not keep_aln(aln, args['mapid'], args['readq'], args['mapq'], args['aln_cov']):
			continue
		mapped += 1
		
		# count alignment
		aln_len = len(aln.query_alignment_sequence)
		gene = db[species_id].genes[gene_id]
		gene.reads += 1
		gene.bases += aln_len
		gene.depth += 1.0*aln_len/gene.length
	
	print("  total aligned reads: %s" % aligned)
	print("  aligned reads after filtering: %s" % mapped)

def quantify_abundance():
	total_depth = 0
	for id, sp in db.items():
		genes = sp.genes.values()
		sp.depth = sum([g.depth for g in genes])
		sp.reads = sum([g.reads for g in genes])
		sp.fract = np.mean([1 if g.reads > 0 else 0 for g in genes])
		total_depth += sp.depth
	for id, sp in db.items():
		sp.abun = 100*sp.depth/total_depth if total_depth > 0 else 0.0
		sp.pres = 1 if sp.fract >= args['presabs_cutoff'] else 0
	total_reads = sum([sp.reads for sp in db.values()])
	print("  mapped reads: %s" % total_reads)
	total_species = sum([sp.pres for sp in db.values()])
	print("  detected species (using presabs_cutoff): %s" % total_species)
	total_species = sum([1 for sp in db.values() if sp.abun > 0])
	print("  detected species (with non-zero abundance): %s" % total_species)

def write_profile():
	out = open('%s/species_profile.tsv' % args['outdir'], 'w')
	fields = ['species_id', 'length', 'genes', 'fract', 'reads', 'depth', 'abund', 'present']
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
		row.append(sp.pres)
		out.write('\t'.join([str(_) for _ in row])+'\n')
	out.close()

def log_time(program_start, module_start):
	current_time = time.time()
	program_time = round(current_time - program_start, 2)
	module_time = round(current_time - module_start, 2)
	peak_ram = round(utility.max_mem_usage(), 2)
	print("  module: %s seconds, program: %s seconds, peak RAM: %s GB" % (module_time, program_time, peak_ram))
	
if __name__ == "__main__":

	import time, utility
	program_start = time.time()
	
	args = parse_arguments()

	check_args(args)
	
	print("\n## Initializing database")
	module_start = time.time()
	db = init_db_info()
	log_time(program_start, module_start)
	
	print("\n## Aligning reads")
	module_start = time.time()
	map_reads_bt2(args)
	log_time(program_start, module_start)

	print("\n## Counting mapped reads")
	module_start = time.time()
	count_reads_bt2(args)
	log_time(program_start, module_start)

	print("\n## Estimating abundance")
	module_start = time.time()
	quantify_abundance()
	log_time(program_start, module_start)

	print("\n## Writing species profile")
	module_start = time.time()
	write_profile()
	log_time(program_start, module_start)




