#!/usr/bin/env python

import os, sys, subprocess as sp, pysam, numpy as np, time, shutil

class Taxon:
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
		description="gMAGsearch: metagenome profiling of species and subspecies from the global gut MAG dataset"
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
	io.add_argument('-d', type=str, dest='db', default=os.environ['MAG_DB'] if 'MAG_DB' in os.environ else None,
		help="""Path to reference database
By default, the MAG_DB environmental variable is used""")
	io.add_argument('--skip', action='store_true', default=False,
		help=argparse.SUPPRESS)
	io.add_argument('-w', dest='weight', choices=['score', 'intra_freq', 'inter_freq', 'none'], default=None,
		help="weight genes based on their score, intra-clade frequency, or inter-clade frequency")
	io.add_argument('--min_intra_freq', type=float, default=0.0)
	io.add_argument('--max_inter_freq', type=float, default=100.0)

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
	""" Check validity of command line arguments """
	# check database
	utility.check_database(args)
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
	
	for line in open('%s' % args['db']):
		
		# parse lines
		if line[0] != '>': continue
		gene_id = line.split()[0].lstrip('>')
		taxon, cluster, genome, gene = gene_id.split('|')
		info = dict([(_.split(':')[0], float(_.split(':')[1])) for _ in line.split(' ', 1)[1].split()])
		
		# get ids
		if taxon.split('_')[0] == 'subspecies':
			spid = taxon.split('_')[1]
			ssid = taxon.split('_', 1)[1]
		else:
			spid = taxon.split('_')[1]
			ssid = None
		
		# init objects
		if spid not in db:
			db[spid] = Taxon()
			db[spid].subspecies = {}
		if (ssid
				and ssid not in db[spid].subspecies):
			db[spid].subspecies[ssid] = Taxon()

		# filter genes
		if info['intra_freq'] < args['min_intra_freq']:
			continue
		elif info['inter_freq'] > args['max_inter_freq']:
			continue

		# store gene
		gene = Gene()
		gene.length = int(info['length'])
		gene.intra_freq = info['intra_freq']
		gene.inter_freq = info['inter_freq']
		gene.score = info['score']
		if ssid:
			db[spid].subspecies[ssid].genes[gene_id] = gene
		else:
			db[spid].genes[gene_id] = gene

	# compute some summary statistics
	for spid, sp in db.items():
		sp.num_genes = len(sp.genes)
		if sp.num_genes > 0:
			sp.max_score = max([g.score for g in sp.genes.values()])
			sp.length = sum([g.length for g in sp.genes.values()])
		for ss in sp.subspecies.values():
			ss.num_genes = len(ss.genes)
			if ss.num_genes > 0:
				ss.max_score = max([g.score for g in ss.genes.values()])
				ss.length = sum([g.length for g in ss.genes.values()])

	# compute gene weights
	for spid, sp in db.items():
		if sp.num_genes > 0:
			weights = fetch_gene_weights(sp.genes.values(), args['weight'])
			for g, w in zip(sp.genes.values(), weights):
				g.weight = w
		for ss in sp.subspecies.values():
			if ss.num_genes > 0:
				weights = fetch_gene_weights(ss.genes.values(), args['weight'])
				for g, w in zip(ss.genes.values(), weights):
					g.weight = w

	return db

def fetch_gene_weights(genes, method):
	if method == 'score':
		values = [max(g.score, 0.0) for g in genes]
	elif method == 'intra_freq':
		values = [g.intra_freq for g in genes]
	elif method == 'inter_freq':
		values = [100.0 - g.inter_freq for g in genes]
	else:
		values = [1.0 for g in genes]
	return [v/sum(values) for v in values] if sum(values) > 0 else [1.0/len(genes) for g in genes]

def map_reads(args):
	import subprocess
	
	# remove existing output
	if os.path.exists('%s/mapped_reads.m8' % args['outdir']):
		os.remove('%s/mapped_reads.m8' % args['outdir'])
	if os.path.exists('%s/mapped_reads.m8.gz' % args['outdir']):
		os.remove('%s/mapped_reads.m8.gz' % args['outdir'])
	
	# build command
	command = 'python %s/stream_seqs.py' % os.path.dirname(os.path.abspath(sys.argv[0]))
	command += ' -1 %s' % args['m1'] # fasta/fastq
	if args['m2']: command += ' -2 %s' % args['m2'] # mate
	if args['max_reads']: command += ' -n %s' % args['max_reads'] # number of reads
	command += ' | hs-blastn align'
	command += ' -query /dev/stdin'
	command += ' -max_target_seqs 1'
	command += ' -num_alignments 1'
	command += ' -db %s' % args['db']
	command += ' -outfmt 6'
	command += ' -num_threads %s' % args['threads']
	command += ' -evalue 1e-3'
	command += ' -out %s/mapped_reads.m8' % args['outdir']

	# run command
	if args['verbose']: print("RUNNING: %s" % command)
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	if process.returncode != 0:
		err_message = "\nError encountered executing:\n%s\n\nError message:\n%s\n" % (command, err)
		sys.exit(err_message)
	if args['verbose'] and len(out) > 0: print("hsblastn stdout: %s" % out)
	if args['verbose'] and len(err) > 0: print("hsblastn stderr: %s" % err)

	# compress output
	if args['verbose']: print("RUNNING: gzip %s/mapped_reads.m8" % args['outdir'])
	process = subprocess.Popen('gzip %s/mapped_reads.m8' % args['outdir'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	if args['verbose'] and len(out) > 0: print("gzip out: %s" % out)
	if args['verbose'] and len(err) > 0: print("gzip err: %s" % err)

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
	import gzip
	formats = [str,str,float,int,float,float,float,float,float,float,float,float]
	fields = ['query','target','pid','aln','mis','gaps','qstart','qend','tstart','tend','evalue','score']
	for line in gzip.open(inpath):
		values = line.rstrip().split()
		yield dict([(field, format(value)) for field, format, value in zip(fields, formats, values)])

def count_mapped_reads(args):
	last_query = None
	alnfile = parse_hsblastn('%s/mapped_reads.m8.gz' % args['outdir'])
	for index, aln in enumerate(alnfile):
		
		# only process each read 1x
		if aln['query'] == last_query:
			continue
		else:
			last_query = aln['query']
	
		# filter alignments
		if not keep_m8_aln(aln, args['mapid'], args['aln_cov']):
			continue

		type, taxon_id = aln['target'].split('|')[0].split('_', 1)
		
		if type == 'subspecies':
			spid = taxon_id.split('_')[0]
			ssid = taxon_id
			taxon = db[spid].subspecies[ssid]
		else:
			spid = taxon_id
			taxon = db[spid]
		
		if aln['target'] not in taxon.genes:
			continue
		else:
			gene = taxon.genes[aln['target']]

		gene.reads += 1
		gene.bases += aln['aln']
		gene.depth += 1.0*aln['aln']/gene.length

def quantify_species():
	total_depth = 0
	for spid, sp in db.items():
		genes = sp.genes.values()
		sp.depth = sum([g.weight * g.depth for g in genes])
		sp.reads = sum([g.reads for g in genes])
		sp.fract = sum([g.weight for g in genes if g.reads > 0])
		total_depth += sp.depth
	for spid, sp in db.items():
		sp.abun = 100*sp.depth/total_depth if total_depth > 0 else 0.0

def write_species_profile():
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

def quantify_subspecies():
	for spid, sp in db.items():
		if len(sp.subspecies) > 0:
			total_depth = 0
			for ss in sp.subspecies.values():
				genes = ss.genes.values()
				ss.depth = sum([g.depth * g.weight for g in genes])
				ss.reads = sum([g.reads for g in genes])
				ss.fract = sum([g.weight for g in genes if g.reads > 0])
				total_depth += ss.depth
			for ss in sp.subspecies.values():
				ss.abun = sp.abun*ss.depth/total_depth if total_depth > 0 else 0.0

def write_subspecies_profile():
	out = open('%s/subspecies_profile.tsv' % args['outdir'], 'w')
	fields = ['taxon_id', 'length', 'genes', 'fract', 'reads', 'depth', 'abund']
	out.write('\t'.join(fields)+'\n')
	for spid, sp in db.items():
		if len(sp.subspecies) > 0:
			for ssid, ss in sp.subspecies.items():
				row = []
				row.append(ssid)
				row.append('%s' % ss.length)
				row.append('%s' % ss.num_genes)
				row.append('%s' % ss.fract)
				row.append('%s' % ss.reads)
				row.append('%s' % ss.depth)
				row.append(ss.abun)
				out.write('\t'.join([str(_) for _ in row])+'\n')
	out.close()

def write_taxon_profile():
	out = open('%s/taxon_profile.tsv' % args['outdir'], 'w')
	fields = ['taxon_id', 'length', 'genes', 'fract', 'reads', 'depth', 'abund']
	out.write('\t'.join(fields)+'\n')
	for spid, sp in db.items():
		if len(sp.subspecies) > 0:
			for ssid, ss in sp.subspecies.items():
				row = []
				row.append(ssid)
				row.append('%s' % ss.length)
				row.append('%s' % ss.num_genes)
				row.append('%s' % ss.fract_covered)
				row.append('%s' % ss.reads)
				row.append('%s' % ss.depth)
				row.append(ss.abun)
				out.write('\t'.join([str(_) for _ in row])+'\n')
		else:
			row = []
			row.append(spid)
			row.append(sp.length)
			row.append(sp.num_genes)
			row.append(sp.fract_covered)
			row.append(sp.reads)
			row.append(sp.depth)
			row.append(sp.abun)
			out.write('\t'.join([str(_) for _ in row])+'\n')
	out.close()

if __name__ == "__main__":

	import time
	start = time.time()
	
	args = parse_arguments()

	import utility
	check_args(args) # --> make sure bowtie2 and/or hsblastn is on PATH

	# initialize database object
	if args['verbose']:
		print("\n## Initializing db")
	db = init_db_info()

	# map reads
	if args['verbose']:
		print("\n## Mapping reads")
	if not (args['skip'] and
			os.path.exists('%s/mapped_reads.m8.gz' % args['outdir'])):
		map_reads(args)

	# count reads
	if args['verbose']:
		print("\n## Counting mapped reads")
	count_mapped_reads(args)

	# summarize mapping
	quantify_species()
	quantify_subspecies()

	# write results
	write_species_profile()
	write_subspecies_profile()

	if args['verbose']:
		print("\n## Pipeline complete")
		print("Time: %s seconds" % round(time.time()-start,2))
		print("RAM: %s GB" % utility.max_mem_usage())




