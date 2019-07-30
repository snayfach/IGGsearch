# IGGsearch

This repository allows you to perform metagenomic species profiling described by http://dx.doi.org/10.1038/s41586-019-1058-x. The pipeline is being actively improved and may change in the future.

<b>Database information</b>

Species-specific marker genes were identified based on two main criterea: i) conserved across genomes within a species, and ii) rarely occurring in genomes from different species. Marker genes for gut species were additionally refined based on co-variation in abundance across metagenomes.

Marker genes were identified from 209,320 genomes from the [IGGdb](https://github.com/snayfach/IGGdb), corresponding to 23,790 species based on 95% genome-wide average nucleotide identity. This genome set includes a large number of genomes from the Human Gut MAG dataset (N=60,664) as well as 148,656 reference genomes and MAGs from PATRIC and IMG.

## Quick start

<b>Install dependencies:</b> 
 
* python (2.7 or 3.6)
	* numpy (1.15.0)
	* pysam (0.15.0) 	
* [Bowtie 2 (2.3.2)](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3)
* [samtools (1.4)](https://github.com/samtools/samtools/releases)

Note: tested versions are indicated in parenthesis, but other versions might also work


Python dependencies can be installed using pip:
`pip install numpy pysam` 

Bowtie2 and samtools should be present on your path:  
`bowtie2 --help`  
`samtools --help`  

<b>Clone the repo from github:</b>
`git clone https://github.com/snayfach/IGGsearch`

<b>Update your environment:</b>   
`export PYTHONPATH=$PYTHONPATH:/path/to/IGGsearch/iggsearch`  
`export PATH=$PATH:/path/to/IGGsearch`  

Note: replace `/path/to` with the correct file path on your system

<b>Download and install the reference database:</b>  

Download: `run_iggsearch.py download --gut-only`  
Unpack: `tar -zxvf iggdb_v1.0.0_gut.tar.gz`   
Install:`export IGG_DB=/path/to/iggdb_v1.0.0_gut`  

Note: use `run_iggsearch.py download` to download full database including non-gut species

<b>Test the code using the dummy dataset:</b>  
`run_iggsearch.py search --m1 test.fastq.gz --outdir test --test`

## Program options

Currently, there are four main modules, which can be viewed using: `run_iggsearch.py -h`

```
Description: Metagenomic species profiling with enhanced coverage of the human gut microbiome

Usage: iggsearch.py <command> [options]

Commands:
   download download reference database of marker genes
     search estimate species abundance from a single metagenome
      merge generate matrix files from multiple runs
   reformat change the file format of from a single run

Note: use iggsearch.py <command> -h to view usage for a specific command
```

Options for the main module can be viewed using: `run_iggsearch.py search -h` 

```
IGGsearch: estimate species abundance from a single metagenome

optional arguments:
  -h, --help            show this help message and exit

input/output:
  --outdir PATH         Directory to store results.
                        Name should correspond to unique identifier for your sample
  --m1 PATH             FASTA/FASTQ file containing 1st mate if using paired-end reads.
                        Otherwise FASTA/FASTQ containing unpaired reads.
                        Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
                        Use comma ',' to separate multiple input files (ex: -1 file1.fq,file2.fq)
  --m2 PATH             FASTA/FASTQ file containing 2nd mate if using paired-end reads.
  --db_dir PATH         Path to reference database. By default, the IGG_DB environmental variable is used

pipeline speed/throughput:
  --max-reads INT       Number of reads to use from input file(s) (use all)
  --threads INT         Number of threads to use for read-alignment (1)
  --no-align            Skip read alignment if 'mapped_reads.bam' already exists (False)
                        Useful for rerunning pipeline with different options
  --test                Perform a quick testing run (False)

read alignment/quality control:
  --mapid FLOAT         Minimum DNA alignment identity between read and marker gene database (95.0)
  --aln_cov FLOAT       Minimum fraction of read covered by alignment (0.75)
  --readq FLOAT         Minimum average base quality score of reads (20.0)
  --mapq FLOAT          Minimum mapping quality of reads (30.0)

species reporting:
  --min-reads-gene INT  Minimum # of reads for detecting marker genes (2)
  --min-perc-genes INT  Minimum % of marker genes detected to report species (40)
  --min-sp-quality INT  Minimum quality score to report species (50)
                        where quality score = completeness - (5 x contamination) of best genome
  --all-species         Presets: --min-reads-gene=0 --min-perc-genes=0 --min-sp-quality=0
  --very-lenient        Presets: --min-reads-gene=1 --min-perc-genes=1 --min-sp-quality=0
  --lenient             Presets: --min-reads-gene=1 --min-perc-genes=15 --min-sp-quality=25
  --strict              Presets: --min-reads-gene=2 --min-perc-genes=40 --min-sp-quality=50 (default)
  --very-strict         Presets: --min-reads-gene=5 --min-perc-genes=60 --min-sp-quality=75
```                        
