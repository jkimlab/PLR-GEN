#!/bin/bash
SRC_PATH=$1 ## script dir path
BIN_PATH=$2 ## binary dir path
ALN=$3  # bam file
REF_FA=$4 # reference sequence fasta file
MAPQ=$5 # minimum mapping quality (cutoff)
PERC=$6  # cutoff percentage for aligned read depth distribution (0 - 100)
MIN_LEN=$7 # minimum length of container
MIN_CNT=$8 # minimum alignment depth for all nodes
BUBBLE_MINOR=$9
OUTDIR=${10} # output directory
PLRID=${11} # psuedo-long read ID

mkdir -p $OUTDIR/raw
echo [$PLRID] Piling-up the alignment
$BIN_PATH/samtools mpileup -q $MAPQ -f $REF_FA --output-QNAME --ff UNMAP,QCFAIL,DUP,SECONDARY $ALN > $OUTDIR/raw/aln.txt
echo [$PLRID] Building graph
perl $SRC_PATH/make_graph.pl $OUTDIR/raw/aln.txt $MIN_LEN $MIN_CNT $OUTDIR/raw
echo [$PLRID] Mapping check 
if [ -f "$OUTDIR/raw/.no_mapping" ]; then
		echo [$PLRID] there is no mapping information - skipped
		rm -rf $OUTDIR
		exit 0
else
		echo [$PLRID] Calculaing bubble stats
		perl $SRC_PATH/Stat_bubble.pl $OUTDIR/raw/graph.out > $OUTDIR/raw/BUBBLE.stat.txt 2> $OUTDIR/raw/BUBBLE.count.txt
		echo [$PLRID] Filtering graph
		perl $SRC_PATH/filtering_bubble_node.pl $OUTDIR/raw/graph.out $OUTDIR/raw/BUBBLE.stat.txt $PERC $BUBBLE_MINOR $OUTDIR > $OUTDIR/processed_graph.out 2> $OUTDIR/log.graph_filtering.txt 
		echo [$PLRID] Merging normal nodes
		$BIN_PATH/bedtools merge -nms -i $OUTDIR/NORMAL.graph.out > $OUTDIR/merge.NORMAL.graph.out
		sed 's/\;//g' $OUTDIR/merge.NORMAL.graph.out > $OUTDIR/merge.NORMAL.nodes.out
		echo [$PLRID] Constructing bubble path
		perl $SRC_PATH/bubble_traverse.pl $OUTDIR/BUBBLE.graph.out > $OUTDIR/consensus_bubble_path.txt 2> $OUTDIR/log.bubble_traverse.txt
		echo [$PLRID] Generating psuedo-long read sequences
		perl $SRC_PATH/longread_seq_generation.pl $OUTDIR/merge.NORMAL.nodes.out $OUTDIR/BUBBLE.graph.out $OUTDIR/consensus_bubble_path.txt $MIN_LEN $OUTDIR 2> $OUTDIR/log.plr_seq.txt
		rm -f $OUTDIR/raw/*.bam $OUTDIR/raw/index.*
		gzip $OUTDIR/raw/*
		echo Done.
fi
