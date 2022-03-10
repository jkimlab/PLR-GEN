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

- Download and install using PLR-GEN package from this github. You can install all dependencies and perl libraries automatically using `build.pl` in this packge. For using this code, `wget`, `unzip`, `make`, `cpanm`, and some perl libraries (FindBin, File::Basename) are needed.

		git clone https://github.com/jkimlab/PLR-GEN.git
		cd PLR-GEN
		./build.pl install
	
- Or, you can prepare all dependencies and perl libraries through creating conda env. 

		git clone https://github.com/jkimlab/PLR-GEN.git
		cd PLR-GEN
		conda env create -f plrgen_env.yml
		conda activate plrgen_env
		./build.pl install

### Docker installation

- If you use Docker, you can download and use PLR-GEN with all dependencies through pulling docker image. 

		docker pull jkimlab/plrgen

### Manual installation

- If you can install PLR-GEN and all its dependencies, you need to prepare all third-party programs into /LOCAL_DISK/PLR-GEN/bin/ directory or make a link of path of installed programs to /LOCAL_DISK/PLR-GEN/bin/ directory. You also need to prepare perl libraries mentioned above (see :point_right: [REQUIREMENTS](https://github.com/jkimlab/PLR-GEN/blob/master/README.md#requirements)). 
 

### (Additional) TAMA installation

- For running PLR-GEN, if you want to use TAMA with `-tama` option for preparation of reference genomes, you can install TAMA using below commands. It will automatically download and install TAMA into the given path, and set ready-made species-level TAMA databases. These TAMA databases require total 300GB disk space, so please carefully set up the path for installation of TAMA. If you do not specify a path for TAMA, the TAMA package and TAMA databases will be set inside "bin" directory in the PLR-GEN. 

1. Installation of TAMA using the code in this package

		/INSTALL_PATH/PLR-GEN/src/TAMA_install.pl PATH_to_INSTALL_for_TAMA
		
2. Installation of TAMA with Docker

		docker run --rm -v /PATH/TO/TAMA_DIR:/tama_dir -t jkimlab/plrgen:latest TAMA_install.pl /tama_dir
		
3. Manual installation

- You can use manually installed TAMA. To use `-tama` option for PLR-GEN, you need to install TAMA and ready-made species-level databases into /LOCAL_DISK/PLR-GEN/bin/ directory or make a link of path of installed TAMA directory to /LOCAL_DISK/PLR-GEN/bin/ directory.

For more information of TAMA, see :point_right: [TAMA github page](https://github.com/jkimlab/TAMA)


## USAGE 
### Input preparation 
1. Input single-end or paired-end NGS short reads
   - Read sequence files are mandatory with fastq format (or compressed by .gz file)
   - You can use single-end reads with `-s` (incompatible with `-1` and `-2`)
   - Or you can use paired-end reads with `-1` and `-2` (incompatible with `-s`)

2. A list of reference genome sequences
   - Reference genome sequence files are mendatory with fasta format (or compressed by .gz file)
   - You can use `-r` or `-ref` option (incompatible with `-tama`)with a list of all the file paths of reference genomes to set up the references for running PLR-GEN

3. (Additional) If TAMA is installed, you can predict putative bacterial species in your dataset using TAMA, running PLR-GEN using `-tama` option (incompatible with `-r`). 

### Running options of PLR-GEN

* Running PLR-GEN

		PLR-GEN.pl [options] -1 <pe1> -2 <pe2> (or -s <se>) -r <ref_list> -o <out_dir>
		

1. Input options (:sparkle: : mandatory option)

   - :sparkle: `-s` : an input sequence file with unpaired reads (single-end reads); a fastq or fastq.gz file
   - :sparkle: `-1` and `-2` : input sequence files with paired-end reads; fastq or fastq.gz files 
   - :sparkle: `-r` or `-ref` : a list with all file paths of reference genome sequences; a text file
   - :sparkle: `-tama` : if this option is used, TAMA program will be used for reference preparation instead of `-r` (default: not used)
   - `-sampling` : if this option is used, reference genomes will be randomly selected and used (default: not used, range: 0-1)

2. Output options 

   - `-o` : a path of output directory (default: ./PR_out; final output: ./PR_out/pseudo-long-read.fa)
   - `-t` : if this option is used, all intermediated output files are saved (default: not used)

3. Running options
   - `-p` : the number of threads
   - `-l` or `-min_length` : minimum length cutoff of generated pseudo-long reads (PLRs) (default: 100)
   - `-q` or `-mapq` : cutoff of mapping quality (used in samtools) (default: 20)
   - `-c` or `-min_count` : cutoff of mapping depth for generating normal nodes and bubbles; after piling-up of read mapping, alignment column less than the cutoff value will be discarded (default: 1)
   - `-d` or `-min_depth` : cutoff of mapping depth of bubbles for filtering bubble nodes; a bubble with less than the cutoff % of mapping depth will be converted to a normal node (default: 1, range: 0-100)

* Help page of PLR-GEN 

		Usage: PLR-GEN.pl [options] -1 <pe1> -2 <pe2> (or -s <se>) -r <ref_list> -o <out_dir>

		== MANDATORY 
		-s	<se>	File with unpaired reads [incompatible with -1 and -2]
		-1	<pe1>	File with #1 mates (paired 1) [incompatible with -s]
		-2	<pe2>	File with #2 mates (paired 2) [incompatible with -s]
		-r|-ref	<ref_list>	The list of reference genome sequence files
		-tama	Reference preparation using TAMA [incompatible with -r|-ref]
		-sampling	<proportion>	proportion to random sampling for references (default: off, range: 0-1)
		-o	<out_dir>	Output directory (default: ./PR.out)

		==Running and filtering options
		-p|-core	<integer>	The number of threads (default: 1)
		-q|-mapq	<integer>	Minimum mapping quality (default: 20)
		-l|-min_length	<integer>	Cutoff of minimum length of pseudo-long reads (default: 100bp)
		-c|-min_count	<integer>	Cutoff of minimum mapping depth for each node (default: 1)
		-d|-min_depth	<integer> Cutoff of mapping depth of bubbles (default: 1, range: 0-100)
				0: all bubbles are used.
				1: bubbles with less than 1% mapping depth from mapping depth distribution of bubbles are converted to normal nodes.
				100: all bubbles are converted to normal nodes

		==Other options
		-t|-temp	If you use -t option, all intermediate files are left.
			Please careful to use this option because it has to be needed very large space.
		-h|-help	Print help page.



## Examples of running PLR-GEN
### Simple examples

1. Using single-end reads and a list of reference genomes

		./PLR-GEN.pl -s s_read.fq -r reference_list.txt -o outdir 
	
2. Using paired-end reads and a list of reference genomes

		./PLR-GEN.pl -1 read_1.fq -2 read_2.fq -r reference_list.txt -o outdir
		
3. Using paired-end reads and a predicted set of references with TAMA

		./PLR-GEN.pl -1 read_1.fq -2 read_2.fq -tama -o outdir
		
4. Using paired-end reads and 10% of reference genomes with random sampling

		./PLR-GEN.pl -1 read_1.fq -2 read_2.fq -r reference_list.txt -sampling 0.1
		
5. Using paired-end reads and 20% of predicted reference genomes with TAMA 

		./PLR-GEN.pl -1 read_1.fq -2 read_2.fq -tama -sampling 0.2

6. Discarding shorter PLRs (< 1000bp)

		./PLR-GEN.pl -1 read_1.fq -2 read_2.fq -tama -min_length 1000
		
7. Using strict mapping cutoff (>= 40)

		./PLR-GEN.pl -1 read_1.fq -2 read_2.fq -tama -mapq 40
		
### Examples of running with the Docker image

- Suppose (1) you pulled the PLR-GEN docker image (jkimlab/plrgen), (2) read_1.fq, read_2.fq are in /LOCAL_DISK/DATA directory, and (3) all reference sequence files are in /LOCAL_DISK/DATA/fasta/ directory. In this situation, you can run PLR-GEN with docker image as followed steps,
	1. Make a reference list file in the /LOCAL_DISK/DATA directory with the file system for the docker image, for an example,
		
		file: /DATA/reference_list.txt
		
			/data_dir/fasta/REF_1.fa
			/data_dir/fasta/REF_2.fa
			/data_dir/fasta/REF_3.fa

	2. Run the docker image, mounting /DATA directory to /data_dir directory of docker image, for an example,

			docker run -v /LOCAL_DISK/DATA:/data_dir -t jkimlab/plrgen PLR-GEN.pl -1 /data_dir/read_1.fq -2 /data_dir/read_2.fq -r /data_dir/reference_list.txt -o /data_dir/output
			
- If you want to use TAMA installed at /LOCAL_DISK/PROGRAMS/TAMA, you need to mount the TAMA directory to /work_dir/bin/ directory of docker image, for an example, 

		docker run -v /LOCAL_DISK/PROGRAMS/TAMA:/work_dir/bin/TAMA -v /LOCAL_DISK/DATA:/data_dir -t jkimlab/plrgen PLR-GEN.pl -1 /data_dir/read_1.fq -2 /data_dir/read_2.fq -r /data_dir/reference_list.txt -o /data_dir/output



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

