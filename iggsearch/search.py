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
		self.alns = []

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
	io.add_argument('--outdir', dest='outdir', type=str, required=True, metavar='PATH',
		help="""Directory to store results.
Name should correspond to unique identifier for your sample""")
	io.add_argument('--m1', type=str, dest='m1', required=True, metavar='PATH',
		help="""FASTA/FASTQ file containing 1st mate if using paired-end reads.
Otherwise FASTA/FASTQ containing unpaired reads.
Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
Use comma ',' to separate multiple input files (ex: -1 file1.fq,file2.fq)""")
	io.add_argument('--m2', type=str, dest='m2', metavar='PATH',
		help="""FASTA/FASTQ file containing 2nd mate if using paired-end reads.""")
	io.add_argument('--db_dir', type=str, default=os.environ['IGG_DB'] if 'IGG_DB' in os.environ else None, metavar='PATH',
		help="""Path to reference database. By default, the IGG_DB environmental variable is used""")

	speed = parser.add_argument_group('pipeline speed/throughput')
	speed.add_argument('--max-reads', type=int, dest='max_reads', metavar='INT',
		help='Number of reads to use from input file(s) (use all)')
	speed.add_argument('--threads', default=1, metavar='INT',
		help='Number of threads to use for read-alignment (1)')
	speed.add_argument('--no-align', action='store_true', default=False,
		help="""Skip read alignment if 'mapped_reads.bam' already exists (False)
Useful for rerunning pipeline with different options""")
	speed.add_argument('--test', action='store_true', default=False,
		help="""Perform a quick testing run (False)""")
	speed.add_argument('--max_genes', type=int, help=argparse.SUPPRESS)

	map = parser.add_argument_group('read alignment/quality control')
	map.add_argument('--mapid', type=float, metavar='FLOAT',
		default=95.0, help='Minimum DNA alignment identity between read and marker gene database (95.0)')
	map.add_argument('--aln_cov', type=float, metavar='FLOAT',
		default=0.75, help='Minimum fraction of read covered by alignment (0.75)')
	map.add_argument('--readq', type=float, metavar='FLOAT',
		default=20.0, help='Minimum average base quality score of reads (20.0)')
	map.add_argument('--mapq', type=float, metavar='FLOAT',
		default=30.0, help='Minimum mapping quality of reads (30.0)')

	species = parser.add_argument_group('species reporting')
	species.add_argument('--min-reads-gene', type=int, default=2, metavar='INT',
		help="""Minimum # of reads for detecting marker genes (2)""")
	species.add_argument('--min-perc-genes', type=int, default=40, metavar='INT',
		help="""Minimum %% of marker genes detected to report species (40)""")
	species.add_argument('--min-sp-quality', type=int, default=50, metavar='INT',
		help="""Minimum quality score to report species (50)
where quality score = completeness - (5 x contamination) of best genome""")
	species.add_argument('--all-species', action='store_true', default=False,
		help="""Presets: --min-reads-gene=0 --min-perc-genes=0 --min-sp-quality=0""")
	species.add_argument('--very-lenient', action='store_true', default=False,
		help="""Presets: --min-reads-gene=1 --min-perc-genes=1 --min-sp-quality=0""")
	species.add_argument('--lenient', action='store_true', default=False,
		help="""Presets: --min-reads-gene=1 --min-perc-genes=15 --min-sp-quality=25""")
	species.add_argument('--strict', action='store_true', default=False,
		help="""Presets: --min-reads-gene=2 --min-perc-genes=40 --min-sp-quality=50 (default)""")
	species.add_argument('--very-strict', action='store_true', default=False,
		help="""Presets: --min-reads-gene=5 --min-perc-genes=60 --min-sp-quality=75""")

	args = vars(parser.parse_args())
	args['file_type'] = iggsearch.utility.auto_detect_file_type(args['m1'].split(',')[0])

	# presets
	if args['all_species']:
		args['min_reads_gene'] = 0
		args['min_perc_genes'] = 0
		args['min_sp_quality'] = 0
	elif args['very_lenient']:
		args['min_reads_gene'] = 1
		args['min_perc_genes'] = 1.0
		args['min_sp_quality'] = 0
	elif args['lenient']:
		args['min_reads_gene'] = 1
		args['min_perc_genes'] = 15.0
		args['min_sp_quality'] = 25
	elif args['strict']:
		args['min_reads_gene'] = 2
		args['min_perc_genes'] = 40.0
		args['min_sp_quality'] = 50
	elif args['very_strict']:
		args['min_reads_gene'] = 5
		args['min_perc_genes'] = 60.0
		args['min_sp_quality'] = 75

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
		error = "\nError: Specified reference database does not exist: %s\n" % args['db_dir']
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
		sp = Species()
		sp.id = r['species_id']
		sp.name = r['species_name']
		sp.gtdb = r['gtdb_taxonomy']
		sp.otus = r['phylo_taxonomy']
		db[r['species_alt_id']] = sp
	
	# add species quality statistics
	inpath = '%s/species_quality.tsv' % os.path.dirname(os.path.realpath(__file__))
	for r in csv.DictReader(open(inpath), delimiter='\t'):
		if r['species_alt_id'] in db:
			sp = db[r['species_alt_id']]
			sp.quality_score = float(r['max_qs'])
	
	# markers
	for index, r in enumerate(csv.DictReader(open('%s.markers' % args['db_base']), delimiter='\t')):
		if args['max_genes'] and index == args['max_genes']: break
		if r['species_alt_id'] not in db: continue
		gene = Gene()
		gene.length = int(r['length'])
		gene.intra_freq = float(r['intra_freq'])
		gene.inter_count = int(r['inter_count'])
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
	if os.path.exists(out) and os.stat(out).st_size > 0 and args['no_align']:
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
	import pysam
	bamfile = pysam.AlignmentFile(bam_path, "r")
	for index, aln in enumerate(bamfile.fetch(until_eof = True)):
	
		# skip species and/or genes that were excluded from the db
		marker_id = bamfile.getrname(aln.reference_id)
		species_id = marker_id.split('|')[0]
		if species_id not in db or marker_id not in db[species_id].genes:
			continue
	
		# filter alignments
		aligned += 1
		if not keep_aln(aln, args['mapid'], args['readq'], args['mapq'], args['aln_cov']):
			continue
		mapped += 1
		
		# count alignment
		aln_len = len(aln.query_alignment_sequence)
		gene = db[species_id].genes[marker_id]
		#gene.alns.append(aln)
		gene.reads += 1
		gene.bases += aln_len
		gene.depth += 1.0*gene.bases/gene.length
	
	print("  total aligned reads: %s" % aligned)
	print("  aligned reads after filtering: %s" % mapped)


def quantify_abundance(args, db):
	
	# summarize alignments per species
	for id, sp in db.items():
		genes = sp.genes.values()
		if sp.length > 0:
			sp.depth = sum([g.bases for g in genes])/float(sp.length)
			sp.reads = sum([g.reads for g in genes])
			sp.detected = len([1 for g in genes if g.reads >= args['min_reads_gene']])
			sp.percent = 100.0*sp.detected/len(genes)

	# prune species
	del_species = []
	for id, sp in db.items():
		if sp.quality_score < args['min_sp_quality']:
			del_species.append(id)
		elif sp.percent < args['min_perc_genes']:
			del_species.append(id)
	for id in del_species:
		del db[id]

	# quantify abundance
	total_depth = 0
	for id, sp in db.items():
		genes = sp.genes.values()
		if sp.length > 0:
			total_depth += sp.depth
	for id, sp in db.items():
		sp.abun = 100*sp.depth/total_depth if total_depth > 0 else 0.0

	print("  reported species: %s" % len(db.values()))
	print("  mapped reads: %s" % sum([sp.reads for sp in db.values()]))


def write_profile(args, db):
	out = open('%s/species_profile.tsv' % args['outdir'], 'w')
	fields = ['species_id', 'species_name',
	          'marker_length', 'marker_count',
			  'percent_markers_detected', 'total_mapped_reads', 'avg_read_depth',
	          'species_abund']
	out.write('\t'.join(fields)+'\n')
	
	species_ids = db.keys()
	abundances = [sp.abun for sp in db.values()]
	species_ids = [id for abun, id in sorted(zip(abundances, species_ids), reverse=True)]
	
	for id in species_ids:
		row = []
		row.append(db[id].id)
		row.append(db[id].name)
		row.append(db[id].length)
		row.append(db[id].num_genes)
		row.append(db[id].percent)
		row.append(db[id].reads)
		row.append(db[id].depth)
		row.append(db[id].abun)
		out.write('\t'.join([str(_) for _ in row])+'\n')
	out.close()

def write_markers(args, db):
	out = open('%s/marker_profile.tsv' % args['outdir'], 'w')
	fields = ['marker_id', 'species_id', 'gene_length', 'intra_freq', 'inter_count', 'num_reads']
	out.write('\t'.join(fields)+'\n')
	for spid, sp in db.items():
		for geneid, gene in sp.genes.items():
			row = [geneid, spid, gene.length, gene.intra_freq, gene.inter_count, gene.reads]
			out.write('\t'.join([str(_) for _ in row])+'\n')

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
	
	print("\n## Identifying species and estimating abundance")
	module_start = time.time()
	quantify_abundance(args, db)
	iggsearch.utility.log_time(program_start, module_start)

	print("\n## Writing species profile")
	module_start = time.time()
	write_profile(args, db)
	iggsearch.utility.log_time(program_start, module_start)



