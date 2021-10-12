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
        
        PLR-GEN.pl [options] -1 <pe1> -2 <pe2> (or -s <se>) -r <ref_list> -o <out_dir>

        == MANDATORY
        -s	<se>	File with unpaired reads
        -1	<pe1>	File with #1 mates (paired 1)
        -2	<pe2>	File with #2 mates (paired 2)
        -r|-ref	<ref_list>	The list of file pathes of reference genome sequences
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
        -h|help	Print help page.

