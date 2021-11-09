PLR-GEN
-----------------

* Generation and application of pseudo-long reads for metagenome assembly


System requirements (tested versions)
-----------------

* Programs

        - perl (v5.22.1)
        - python (3)
        - git (2.7.4)
        - gcc, g++ (7.5.0)
        - make (GNU Make 4.1)
        - zip (3.0)
        - wget (1.17.1)
        - glibc (2.14+) 

* Perl libraries

        - Parallel::ForkManager 
        - Getopt::Long
        - File::Basename
        - Scalar::Util
        - FindBin
        - Math::Round
        
        

Download and installation
-----------------

* Download and install (source code)

        git clone https://github.com/jkimlab/PLR-GEN.git
        cd PLR-GEN
        ./build.pl install


* Download Docker image (with all dependencies)  RECOMMENDED! 

        docker pull mksim/plrgen:latest
        

Running PLR-GEN
-----------------

* Running options
        
        PLR-GEN.pl [options] -1 <pe1> -2 <pe2> (or -s <se>) -r <ref_list> -o <out_dir>

        == MANDATORY
        -s	<se>	file with unpaired reads
        -1	<pe1>	file with #1 mates (paired 1)
        -2	<pe2>	file with #2 mates (paired 2)
        -r|-ref	<ref_list>	The list of file paths of reference genome sequences
        -o	<out_dir>	Output directory (default: ./PLR.out)

        ==Running and filtering options
        -p|-core	<integer>	the number of threads (default: 1)
        -q|-mapq	<integer>	minimum mapping quality (default: 20)
        -l|-min_length	<integer>	cutoff of minimum length of pseudo-long reads (default: 100bp)
        -c|-min_count	<integer>	cutoff of minimum mapping depth for each node (default: 1)
        -d|-min_depth	<integer>       cutoff of mapping depth of bubbles (default: 1, 0-100)
        		0: all bubble are used.
        		1: bubbles with less than lowest 1% mapping depth from distribution of mapping depth of bubbles are converted to normal node.
        		100: all bubbles are converted to normal nodes        
        
        ==Other options
        -t|-temp	If you use -t option, all intermediate files are left.
        	        Please use this option carefully, because it could to be needed very large space.
        -h|-help	Print help page.


* Parameter information

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
       

Example data
-----------------
        
* Input short reads (Human Microbiome Project dataset)
  
  Metagenome short paired-end reads of HMP dataset used to application of PLR-GEN (in manuscript) are available from NCBI SRA.
        
        - Accession number : SRR2822457
        - NCBI SRA URL : https://www.ncbi.nlm.nih.gov/sra/SRR2822457
        
* Reference genome sequences 

  All used reference genome sequence files for running PLR-GEN (in manuscript) are in ref_fa directory of this package.
  There is list.txt file in the ref_fa directoy that is an example of the list of reference file.
        
