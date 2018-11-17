# IGGsearch
<b>Metagenomic species profiling with enhanced coverage of the human gut microbiome</b>

IGGsearch accurately quantifies species presence-absence and species abundance by mapping reads to a database of species-specific marker genes. 

Marker genes were identified based on two main criterea: i) conserved across genomes within a species, and ii) rarely occurring in genomes from different species. Additionally, marker genes for gut species were refined based on co-variation in abundance across metagenomes, which was especially important for species with just a single genome where gene conservation was challenging to estimate.

Marker genes were identified from 209,320 genomes from the [IGGdb](https://github.com/snayfach/IGGdb), corresponding to 23,790 species based on 95% genome-wide average nucleotide identity. This genome set includes a large number of genomes from the Human Gut MAG dataset (N=60,664) as well as 148,656 reference genomes and MAGs from PATRIC and IMG.

Currently, there are two available marker-gene databases: for all 23,790 species and for only 4,558 species from the human gut microbiome. 

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
`export PYTHONPATH=$PYTHONPATH:/path/to/IGGsearch`  
`export PATH=$PATH:/path/to/IGGsearch/scripts`  

Note: replace `/path/to` with the correct file path on your system

<b>Download and install the reference database:</b>  

Download: `run_iggsearch.py download --gut-only`  
Unpack: `tar -zxvf iggdb_v1.0.0_gut.tar.gz`   
Install:`export IGG_DB=/path/to/iggdb_v1.0.0_gut`  

Note: use `run_iggsearch.py download` to download full database including non-gut species

<b>Test the code using the dummy dataset:</b>  
`run_iggsearch.py search -1 test.fastq.gz -o test --test`

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
  -h, --help         show this help message and exit

input/output:
  --outdir PATH      Directory to store results.
                     Name should correspond to unique identifier for your sample
  --m1 PATH          FASTA/FASTQ file containing 1st mate if using paired-end reads.
                     Otherwise FASTA/FASTQ containing unpaired reads.
                     Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
                     Use comma ',' to separate multiple input files (ex: -1 file1.fq,file2.fq
  --m2 PATH          FASTA/FASTQ file containing 2nd mate if using paired-end reads.
  --db_dir PATH      Path to reference database. By default, the IGG_DB environmental variable is used

pipeline speed:
  --max-reads INT    Number of reads to use from input file(s) (use all)
  --threads INT      Number of threads to use for read-aignment (1)
  --no-align         Skip read alignment if 'mapped_reads.bam' already exists (False)
                     Useful for rerunning pipeline with different options
  --test             Perform a quick testing run (False)

alignment/quality control:
  --mapid FLOAT      Discard reads with alignment identity < MAPID (95.0)
  --aln_cov FLOAT    Discard reads with alignment coverage < ALN_COV (0.75)
  --readq FLOAT      Minimum average-base-quality per read (20.0)
  --mapq FLOAT       Minimum map quality score per read (0)

species filtering:
  --all              Output results for all species, including those that were not detected (False)
  --no-sort          Do not order species by abundance in output file (False)
                     Useful when combined with '--all' to enforce same ordering of species across multiple output files
  --hq-only          Only report species with at least 1 high-quality genome (False)
  --min-markers INT  Exclude species with fewer than <min-markers> (0)
  --pres FLOAT       Threshold for determining species presence-absence,
                     defined at the percent of a species' marker genes with >=1 mapped read.
                     Useful for eliminating spurious hits (15)
```                        