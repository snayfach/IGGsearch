# IGGsearch
Quantification of gut lineages from IGGdb from metagenomes

### Installation

<b>Clone the repo from github:</b>
`git clone https://github.com/snayfach/IGGsearch`

<b>Download the reference database:</b>  

From your browser:   
[http://lighthouse.ucsf.edu/IGGdb/iggsearch_db_v1.0.tar.gz](http://lighthouse.ucsf.edu/IGGdb/iggsearch_db_v1.0.tar.gz)

Or, from the command line:   
on Unix: `wget http://lighthouse.ucsf.edu/IGGdb/iggsearch_db_v1.0.tar.gz`  
on OSX: `curl http://lighthouse.ucsf.edu/IGGdb/iggsearch_db_v1.0.tar.gz > iggsearch_db_v1.0.tar.gz`

Unpack the database: `tar -zxvf iggsearch_db_v1.0.tar.gz`   

<b>Install the hs-blastn alignment tool:</b>  
1. `git clone https://github.com/chenying2016/queries.git`  
2. `cd queries/hs-blastn-src`  
3. `make`   
4. `export PATH=$PATH:$PWD` or place the `hs-blastn` binary in a location located on your $PATH


### Testing

To list options, use: `python iggsearch.py -h`  

Run IGGsearch on a very small test dataset
`python iggsearch.py -1 test.fastq.gz -d iggsearch_db_v1.0 -o test --verbose`
