#!/usr/bin/env python

class Species:
	def __init__(self):
		self.genes = {}
		self.reads = 0
		self.bases = 0.0
		self.depth = 0.0
		self.num_genes = 0
		self.num_covered = 0
		self.percent = 0.0
		self.length = 0

class Gene:
	def __init__(self):
		self.length = None
		self.reads = 0
		self.bases = 0
		self.depth = 0.0
		self.intra_freq = None
		self.inter_freq = None

def import_libraries():
	import importlib
	modnames = ["os", "sys", "subprocess", "numpy", "time", "pysam", "csv", "time", "argparse", "iggsearch", "iggsearch.utility"]
	for lib in modnames:
		errors = []
		try:
			globals()[lib] = importlib.import_module(lib)
		except:
			errors.append(lib)
	if len(errors) > 0:
		sys.exit("\nCould not import the following libraries: %s \n" % str(errors))

def parse_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="IGGsearch: estimate species abundance from a single metagenome"
	)
	parser.add_argument('program', help=argparse.SUPPRESS)
	io = parser.add_argument_group('input/output')
	io.add_argument('-o', dest='outdir', type=str, required=True,
		help="""Path to directory to store results.
Directory name should correspond to sample identifier""")
	io.add_argument('-1', type=str, dest='m1', required=True,
		help="""FASTA/FASTQ file containing 1st mate if using paired-end reads.
Otherwise FASTA/FASTQ containing unpaired reads.
Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
Use comma ',' to separate multiple input files (ex: -1 file1.fq,file2.fq""")
	io.add_argument('-2', type=str, dest='m2',
		help="""FASTA/FASTQ file containing 2nd mate if using paired-end reads.""")
	io.add_argument('-d', type=str, dest='db_dir', default=os.environ['IGG_DB'] if 'IGG_DB' in os.environ else None,
		help="""Path to reference database. By default, the IGG_DB environmental variable is used""")
	io.add_argument('--all', action='store_true', default=False,
		help="""Output results for all species, including those that were not detected (False)""")
	io.add_argument('--no-sort', action='store_true', default=False,
		help="""Do not order species by descreasing abundnce in output file (False)
Useful with combined with '--all' to enforce same ordering of species across multiple output files""")
	io.add_argument('--hq-only', action='store_true', default=False,
		help="""Only report results for species with at least 1 high-quality genome (False)""")
	speed = parser.add_argument_group('pipeline speed')
	speed.add_argument('-n', type=int, dest='max_reads',
		help='# reads to use from input file(s) (use all)')
	speed.add_argument('-t', dest='threads', default=1,
		help='Number of threads to use (1)')
	speed.add_argument('--no-align', action='store_true', default=False,
		help="""Skip read alignment if <outdir>/mapped_reads.bam already exists (False)""")
	speed.add_argument('--test', action='store_true', default=False,
		help="""Perform a quick testing run (False)""")
	speed.add_argument('--max_genes', type=int, help=argparse.SUPPRESS)
	map = parser.add_argument_group('alignment/quality control')
	map.add_argument('--mapid', type=float, metavar='FLOAT',
		default=95.0, help='Discard reads with alignment identity < MAPID (95.0)')
	map.add_argument('--aln_cov', type=float, metavar='FLOAT',
		default=0.75, help='Discard reads with alignment coverage < ALN_COV (0.75)')
	map.add_argument('--readq', type=float, metavar='FLOAT',
		default=20.0, help='Minimum average-base-quality per read (20.0)')
	map.add_argument('--mapq', type=float, metavar='FLOAT',
		default=0, help='Minimum map quality score per read (0)')
	map.add_argument('--pres', type=float, default=15,
		help="""Species presence-absence threshold,
defined at the percent of a species' marker genes with >=1 mapped read.
Useful for eliminating spurious hits (15)""")
	map.add_argument('--min-markers', type=int, default=None,
		help="""Exclude species with fewer than <min-markers> (0)""")
	args = vars(parser.parse_args())
	args['file_type'] = iggsearch.utility.auto_detect_file_type(args['m1'].split(',')[0])
	check_args(args)
	return args

def check_args(args):
	# check database
	if args['db_dir'] is None:
		error = "\nError: No reference database specified\n"
		error += "Use the flag -d to specify a database,\n"
		error += "Or set the IGG_DB environmental variable: export IGG_DB=/path/to/igg_db\n"
		sys.exit(error)
	args['db_base'] = '%s/%s' % (args['db_dir'], args['db_dir'].split('/')[-1])
	if not os.path.exists(args['db_dir']):
		error = "\nError: Specified reference database does not exist: %s" % os.path.dirname(args['db_dir'])
		sys.exit(error)
	for file in []:
		path = '%s/%s' % (args['db_dir'], file)
		if not os.path.exists(path):
			error = "\nError: Could not locate required database file: %s\n" % path
			sys.exit(error)
	# check executables
	for exe in ['bowtie2', 'samtools']:
		if not iggsearch.utility.which(exe):
			sys.exit("\nError: required program '%s' not executable or not found on $PATH\n" % exe)
	if args['test']:
		args['max_reads'] = 1000
		args['max_genes'] = 100000
		args['no_align'] = True
		sys.stderr.write("\nWarning: this is a quick test run; results should not be interpreted\n")
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

def init_db_info(args):
	
	db = {}

	# initialize species
	for r in csv.DictReader(open('%s.species' % args['db_base']), delimiter='\t'):
		if args['hq_only'] and r['is_high_quality']=='No':
			continue
		if args['min_markers'] and int(r['marker_count']) < args['min_markers']:
			continue
		sp = Species()
		sp.id = r['species_id']
		sp.name = r['species_name']
		sp.gtdb = r['gtdb_taxonomy']
		sp.otus = r['phylo_taxonomy']
		db[r['species_alt_id']] = sp
	
	# markers
	for index, r in enumerate(csv.DictReader(open('%s.markers' % args['db_base']), delimiter='\t')):
		if args['max_genes'] and index == args['max_genes']: break
		if r['species_alt_id'] not in db: continue
		gene = Gene()
		gene.length = int(r['length'])
		db[r['species_alt_id']].genes[r['marker_id']] = gene

	# compute some summary statistics
	for id, sp in db.items():
		sp.num_genes = len(sp.genes)
		if sp.num_genes > 0:
			sp.length = sum([g.length for g in sp.genes.values()])

	print("  total species: %s" % len(db))
	print("  total genes: %s" % sum([len(sp.genes) for sp in db.values()]))
	
	return db


def map_reads_bt2(args):
	
	out = '%s/mapped_reads.bam' % args['outdir']
	if os.path.exists(out) and args['no_align']:
		print("  nothing to do")
		return
	
	# Run bowtie2
	command = 'bowtie2 --no-unal '
	command += '-f ' if args['file_type'] == 'fasta' else '-q '
	command += '-x %s ' % args['db_base']
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
	out, err = [str(_) for _ in process.communicate()]
	if process.returncode != 0:
		err_message = "\nError encountered executing:\n%s\n\nError message:\n%s\n" % (command, err)
		sys.exit(err_message)

	# write log
	with open('%s/read_alignment.log' % args['outdir'], 'w') as f:
		f.write(out+'\n'+err)

def keep_aln(aln, min_pid, min_readq, min_mapq, min_aln_cov):
	align_len = len(aln.query_alignment_sequence)
	query_len = aln.query_length
	# min pid
	if 100*(align_len-dict(aln.tags)['NM'])/float(align_len) < min_pid:
		return False
	# min read quality
	elif numpy.mean(aln.query_qualities) < min_readq:
		return False
	# min map quality
	elif aln.mapping_quality < min_mapq:
		return False
	# min aln cov
	elif align_len/float(query_len)  < min_aln_cov:
		return False
	else:
		return True

def count_reads_bt2(args, db):
	aligned = 0
	mapped = 0
	bam_path =  '%s/mapped_reads.bam' % args['outdir']
	bamfile = pysam.AlignmentFile(bam_path, "r")
	for index, aln in enumerate(bamfile.fetch(until_eof = True)):
	
		# skip species and/or genes that were excluded from the db
		marker_id = bamfile.getrname(aln.reference_id)
		species_id = marker_id.split('|')[0]
		if species_id not in db or marker_id not in db[species_id].genes:
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
		gene = db[species_id].genes[marker_id]
		gene.reads += 1
		gene.bases += aln_len
		gene.depth += 1.0*gene.bases/gene.length
	
	print("  total aligned reads: %s" % aligned)
	print("  aligned reads after filtering: %s" % mapped)

def quantify_abundance(args, db):
	total_depth = 0
	for id, sp in db.items():
		genes = sp.genes.values()
		if sp.length > 0:
			sp.depth = sum([g.bases for g in genes])/float(sp.length)
			sp.reads = sum([g.reads for g in genes])
			sp.percent = 100*numpy.mean([1 if g.reads > 0 else 0 for g in genes])
			total_depth += sp.depth
	for id, sp in db.items():
		sp.pres = 1 if sp.percent >= args['pres'] else 0
	for id, sp in db.items():
		sp.abun = 100*sp.depth/total_depth if total_depth > 0 else 0.0 #and sp.pres = 1
	total_reads = sum([sp.reads for sp in db.values()])
	print("  mapped reads: %s" % total_reads)
	total_species = sum([sp.pres for sp in db.values()])
	print("  detected species (using presabs_cutoff): %s" % total_species)
	total_species = sum([1 for sp in db.values() if sp.abun > 0])
	print("  detected species (with non-zero abundance): %s" % total_species)

def write_profile(args, db):
	out = open('%s/species_profile.tsv' % args['outdir'], 'w')
	fields = ['species_id', 'species_name',
	          'marker_length', 'marker_count',
			  'percent_markers_detected', 'total_mapped_reads', 'avg_read_depth',
	          'species_abund', 'species_presence']
	out.write('\t'.join(fields)+'\n')
	
	species_ids = db.keys()
	if not args['no_sort']:
		abundances = [sp.abun for sp in db.values()]
		species_ids = [id for abun, id in sorted(zip(abundances, species_ids), reverse=True)]
	
	for id in species_ids:
		if not args['all'] and db[id].reads == 0:
			continue
		row = []
		row.append(db[id].id)
		row.append(db[id].name)
		row.append(db[id].length)
		row.append(db[id].num_genes)
		row.append(db[id].percent)
		row.append(db[id].reads)
		row.append(db[id].depth)
		row.append(db[id].abun)
		row.append(db[id].pres)
		out.write('\t'.join([str(_) for _ in row])+'\n')
	out.close()

def main():

	import_libraries()
	program_start = time.time()
	args = parse_arguments()
	
	print("\n## Initializing database")
	module_start = time.time()
	db = init_db_info(args)
	iggsearch.utility.log_time(program_start, module_start)
	
	print("\n## Aligning reads")
	module_start = time.time()
	map_reads_bt2(args)
	iggsearch.utility.log_time(program_start, module_start)

	print("\n## Counting mapped reads")
	module_start = time.time()
	count_reads_bt2(args, db)
	iggsearch.utility.log_time(program_start, module_start)

	print("\n## Estimating abundance")
	module_start = time.time()
	quantify_abundance(args, db)
	iggsearch.utility.log_time(program_start, module_start)

	print("\n## Writing species profile")
	module_start = time.time()
	write_profile(args, db)
	iggsearch.utility.log_time(program_start, module_start)




