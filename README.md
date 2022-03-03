# Overview

PLR-GEN is a tool for the generation of pseudo-long-reads (PLRs) by using short-reads of a metagenomic sample and microbial reference genome sequences as input. 

## REQUIREMENTS
### Third party programs

- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [BEDtools](https://bedtools.readthedocs.io/en/latest/)
- [SAMtools](http://www.htslib.org/)

### Perl libraries

- Parallel::ForkManager 
- Getopt::Long
- File::Basename
- Scalar::Util
- FindBin
- Math::Round

## INSTALLATION
### Installation using PLR-GEN package

- Download and install using PLR-GEN package from this github. You can install all dependencies and perl libraries automatically using `build.pl` in this packge. For using this code, `wget`, `unzip`, `make`, `cpanm`, and some perl libraries (FindBin, File::Basename)

		git clone https://github.com/jkimlab/PLR-GEN.git
		cd PLR-GEN
		./build.pl install
		./PLR-GEN.pl
	
- Or, you can prepare all dependencies and perl libraries through creating conda env. 

		git clone https://github.com/jkimlab/PLR-GEN.git
		cd PLR-GEN
		conda env create -f plrgen_env.yml
		conda activate plrgen_env

### Docker installation

- If you use Docker, you can download and use PLR-GEN with all dependencies through pulling docker image. 

		docker pull jkimlab/plrgen

### Manual installation

- If you can install PLR-GEN and all its dependencies, you need to prepare all third-party programs and adding to $PATH, also need to prepare perl libraries (see :point_right: [REQUIREMENTS](https://github.com/jkimlab/PLR-GEN/blob/master/README.md#requirements)). 
 

### (Additional) TAMA installation

- If you want to use TAMA for preparation of reference genomes, you can install TAMA using below commands. It automatically download and install TAMA into the given path, and set the ready-made species-level TAMA databases. These TAMA databases require total 300GB disk space, so please carefully set up the path for installation of TAMA. If you do not specify a path for TAMA, the TAMA package and TAMA databases will be set inside "bin" directory in the PLR-GEN.

1. Installation of TAMA using the code in this package

		/LOCAL_PATH/PLR-GEN/src/TAMA_install.pl PATH_to_INSTALL_for_TAMA
		
2. Installation of TAMA with Docker

		docker run --rm -v /PATH/TO/TAMA_DIR:/tama_dir -t jkimlab/plrgen:latest /work_dir/src/TAMA_install.pl /tama_dir

For more information of TAMA, see :point_right: [TAMA github page](https://github.com/jkimlab/TAMA)


## USAGE 
### Input preparation 
1. Input single-end or paired-end NGS short reads
   - Read sequences are mendatory with fastq format
   - You can use single-end reads with `-s` (incompatible with paired-end option)
   - Or you can use paired-end reads with `-1` and `-2` (incompatible with single-end option)

2. A list of reference genome sequences
   - 
### Running options of PLR-GEN
        
	Usage: PLR-GEN.pl [options] -1 <pe1> -2 <pe2> (or -s <se>) -r <ref_list> -o <out_dir>

	== MANDATORY
	-s	<se>	File with unpaired reads [incompatible with -1 and -2]
	-1	<pe1>	File with #1 mates (paired 1) [incompatible with -s]
	-2	<pe2>	File with #2 mates (paired 2) [incompatible with -s]
	-r|-ref	<ref_list>	The list of reference genome sequence files
	-tama	Reference preparation using TAMA [incompatible with -r|-ref]
	-sampling	<proportion>	proportion to random sampling for references (default: off, range: 0-1)
	-o	<out_dir>	Output directory (default: ./PR.out

	==Running and filtering options
	-p|-core	<integer>	the number of threads (default: 1)
	-q|-mapq	<integer>	minimum mapping quality (default: 20)
	-l|-min_length	<integer>	cutoff of minimum length of pseudo-long reads (default: 100bp)
	-c|-min_count	<integer>	cutoff of minimum mapping depth for each node (default: 1)
	-d|-min_depth	<integer> cutoff of mapping depth of bubbles (default: 1, range: 0-100)
			0: all bubbles are used.
			1: bubbles with less than 1% mapping depth from mapping depth distribution of bubbles are converted to normal nodes.
			100: all bubbles are converted to normal nodes

	==Other options
	-t|-temp	If you use -t option, all intermediate files are left.
		Please careful to use this option because it has to be needed very large space.
	-h|-help	Print help page.


### Parameter information

	[ Input sequence files ]

	-s      
		an input sequence file with unpaired reads (single-end reads)
		a fastq or fastq.gz file

	-1, -2 
		input sequence files with paired-end reads
		fastq or fastq.gz files 

	-r|-ref 
		a list with file paths of reference genome sequences
		a text file

	-o 
		a name (or specific path) of output directory
		all outputs will be generated in this directory 

	[ PLR-GEN running options ]

	-p 
		the number of threads
		default: 1

	-q|-mapq
		cutoff of mapping quality (used in samtools)
		read mapping from bowtie outputs less than this cutoff value will be discarded
		default: 20

	-l|-min_length
		minimum length cutoff of generated pseudo-long reads (PLRs)
		generated PLRs less than the cutoff value will be discarded
		default: 100

	-c|-min_count
		cutoff of mapping depth for generating normal nodes and bubbles
		after piling-up of read mapping, alignment column less than the cutoff value will be discarded
		default: 1

	-d|-min_depth
		cutoff of mapping depth of bubbles for filtering bubble nodes
		a bubble with less than the cutoff % of mapping depth will be converted to a normal node

	[ Other options ]
	-t|-temp
		all intermediate output and intermediate log files are left
		It could be needed very large space
        
* Input and output file format

        * a list of reference sequence file 
        
        - absolute or relative paths of reference sequence files
        - example
                [local_data_path]/ref_1.fa
                [local_data_path]/ref_2.fa
                ... 

        * pseudo-long_read.fa
                generated PLR sequence fasta file (final output)
                
        * log.PSEUDO-LONG_READ_GEN.txt
                recorded log messages and used parameter
       

Running PLR-GEN with HMP dataset
-----------------
        
* Input short reads (Human Microbiome Project dataset)
  
  Metagenome short paired-end reads of HMP dataset used to application of PLR-GEN (in manuscript) are available from NCBI SRA.
        
        - Accession number : SRR2822457
        - NCBI SRA URL : https://www.ncbi.nlm.nih.gov/sra/SRR2822457

  If you have SRA toolkit, you can download this dataset using Download.sh file in 'example' directory, using the command below.
  
        cd PLR-GEN/example
        ./Download.sh 
        
        
* Reference genome sequences 

  All used reference genome sequence files for running PLR-GEN (in manuscript) are in ref_fa directory of this package.
  There is list.txt file in 'example' directoy that is an example of the list of reference file.

* Running PLR-GEN with example dataset

  In 'example' directory, there is a file with command lines (named CMD.sh) for running PLR-GEN using the above short paired-end reads and reference genome sequences.
  You need to prepare SRA toolkit program if you want to download paired-end reads using the Download.sh file.
  Or you can download the paired-end reads personally and then locate them in 'example' directory.
  When all files of example dataset (metagenome short reads and a file list of reference genome sequences), use the command below.
  
        cd PLR-GEN/example
        ./CMD.sh <threads> <outdir>

