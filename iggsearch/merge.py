#!/usr/bin/env python


def import_libraries():
	import importlib
	modnames = ["os", "sys", "time", "csv", "time", "argparse", "iggsearch", "iggsearch.utility"]
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
		description="IGGsearch: merge: create matrix files from multiple runs"
	)
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('--outdir', type=str, required=True, metavar='PATH',
		help="""Directory to store results.
Name should correspond to unique identifier for your sample""")
	parser.add_argument('--input', type=str, required=True, metavar='PATH',
		help="""FASTA/FASTQ file containing 1st mate if using paired-end reads.
Otherwise FASTA/FASTQ containing unpaired reads.
Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
Use comma ',' to separate multiple input files (ex: -1 file1.fq,file2.fq""")
	parser.add_argument('--intype', choices=['list','file','dir'], required=True, metavar="CHOICE",
		help="""Choose from 'list', 'dir', 'file'.
list: --intype is a comma-separated list (ex: /samples/sample_1,/samples/sample_2)
dir: --intype is a directory containing all samples (ex: /samples)
file: --intype is a file containing paths to samples (ex: sample_paths.txt)""")
	parser.add_argument('--db_dir', type=str, default=os.environ['IGG_DB'] if 'IGG_DB' in os.environ else None, metavar='PATH',
		help="""Path to reference database. By default, the IGG_DB environmental variable is used""")
	parser.add_argument('--max-samples', type=int,
		help="""Maximum samples to process from --input; useful for quick tests""")

	args = vars(parser.parse_args())
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
	# create output directory
	if not os.path.isdir(args['outdir']):
		os.makedirs(args['outdir'])

def list_input_files(args):
	samples = []
	if args['intype'] == 'list':
		samples = args['input'].split(',')
	elif args['intype'] == 'file':
		if not os.path.isfile(args['input']):
			sys.exit("\nError: specified file does not exist: %s\n" % args['input'])
		samples = [_.rstrip() for _ in open(args['input'])]
	elif args['intype'] == 'dir':
		if not os.path.isdir(args['input']):
			sys.exit("\nError: specified directory does not exist: %s\n" % args['input'])
		samples = ['%s/%s' % (args['input'], _) for _ in os.listdir(args['input'])]
	inpaths = []
	for sample in samples:
		inpath = '%s/species_profile.tsv' % sample
		if os.path.exists(inpath):
			inpaths.append(inpath)
		else:
			sys.stderr.write("  warning: file does not exist: %s" % inpath)
	if args['max_samples']:
		return inpaths[0:args['max_samples']]
	return inpaths


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
	
	return db

def main():

	import_libraries()
	program_start = time.time()
	args = parse_arguments()
	
	print("\n## Listing input")
	module_start = time.time()
	inpaths = list_input_files(args)
	print("  %s total samples" % len(inpaths))
	iggsearch.utility.log_time(program_start, module_start)
	
	print("\n## Reading data")
	module_start = time.time()
	data = {}
	species_ids = set()
	for inpath in inpaths:
		sample = inpath.split('/')[-2]
		data[sample] = {}
		with open(inpath) as f:
			for r in csv.DictReader(f, delimiter='\t'):
				data[sample][r['species_id']] = r
				species_ids.add(r['species_id'])
	iggsearch.utility.log_time(program_start, module_start)
	species_ids = sorted(list(species_ids))

	print("\n## Writing matrix files")
	module_start = time.time()
	fields = [('species_presence', '0'),
	          ('avg_read_depth', '0.0'),
	          ('species_abund', '0.0'),
	          ('percent_markers_detected', '0.0'),
	          ('total_mapped_reads', '0')]
	header = 'sample_id\t'+'\t'.join(species_ids)+'\n'
	for field, default in fields:
		path = '%s/%s.tsv' % (args['outdir'], field)
		with open(path, 'w') as file:
			file.write(header)
			for sample in data.keys():
				file.write(sample)
				for species in species_ids:
					if species in data[sample]:
						value = str(data[sample][species][field])
					else:
						value = default
					file.write('\t'+value)
				file.write('\n')
	iggsearch.utility.log_time(program_start, module_start)
	





