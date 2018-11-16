# IGGsearch
Metagenomic species profiling with enhanced coverage of the human gut microbiome

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


