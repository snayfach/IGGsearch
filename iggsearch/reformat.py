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
		description="IGGsearch: reformat species-level matrices "
	)
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('--outdir', type=str, required=True, metavar='PATH',
		help="""Directory to store results""")
	parser.add_argument('--indir', type=str, required=True, metavar='PATH',
		help="""Input directory containing matrix files""")
	parser.add_argument('--taxdb', metavar='CHOICE', type=str, choices=['otus', 'gtdb'], required=True,
help="""Taxonomy database: 'otus', 'gtdb'
otus - operational taxonomic units based on phylogeny
gtdb - Genome Taxonomy Database
""")
	parser.add_argument('--taxrank', metavar='CHOICE', type=str, choices=['genus', 'family', 'order', 'class', 'phylum'], required=True, help="Taxonomic rank: 'genus', 'family', 'order', 'class', 'phylum'")
	
	parser.add_argument('--db_dir', type=str, default=os.environ['IGG_DB'] if 'IGG_DB' in os.environ else None, metavar='PATH',
		help="""Path to reference database.
By default, the IGG_DB environmental variable is used""")
	parser.add_argument('--max-samples', type=int, metavar='INT',
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

def fetch_species_to_clade(args):
	species_to_clade = {}
	rank_to_index = {'phylum':1, 'class':2, 'order':3, 'family':4, 'genus':5}
	index = rank_to_index[args['taxrank']]
	for r in csv.DictReader(open('%s.species' % args['db_base']), delimiter='\t'):
		if args['taxdb'] == 'gtdb':
			taxonomy = r['gtdb_taxonomy']
		else:
			taxonomy = r['phylo_taxonomy']
		clades = taxonomy.split(';')
		if len(clades) <= index:
			species_to_clade[r['species_id']] = args['taxrank'][0]+'_unclassified'
		else:
			species_to_clade[r['species_id']] = clades[index]
	return species_to_clade

def main():

	import_libraries()
	program_start = time.time()
	args = parse_arguments()
	
	species_to_clade = fetch_species_to_clade(args)
	
	fields = [('avg_read_depth', float, sum),
	         ('species_presence', int, sum),
	         ('percent_markers_detected', float, max),
	         ('species_abund', float, sum),
			 ('total_mapped_reads', int, sum)]
	
	for field, format, fun in fields:
		data = {}
		inpath = '%s/%s.tsv' % (args['indir'], field)
		print ("## Aggregating data for '%s'" % inpath)
		module_start = time.time()
		with open(inpath) as f:
			species_ids = next(f).split()[1:]
			clades = [species_to_clade[id] for id in species_ids]
			for i, l in enumerate(f):
				if args['max_samples'] and i == args['max_samples']:
					break
				sample_id = l.split()[0]
				values = l.split()[1:]
				if sample_id not in data:
					data[sample_id] = {}
				for clade, value in zip(clades, values):
					if clade not in data[sample_id]:
						data[sample_id][clade] = []
					data[sample_id][clade].append(format(value))

		outpath = '%s/%s.tsv' % (args['outdir'], field)
		with open(outpath, 'w') as f:
			clades = sorted(data.values()[0].keys())
			header = 'sample_id\t'+'\t'.join(clades)+'\n'
			f.write(header)
			for sample in data.keys():
				f.write(sample)
				for clade in clades:
					value = fun(data[sample][clade])
					f.write('\t'+str(value))
				f.write('\n')
				
		iggsearch.utility.log_time(program_start, module_start)

