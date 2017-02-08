#!/bin/bash
#$ -l h_cpu=72:00:00
#$ -l h_vmem=8G
#$ -cwd
#$ -N bed_stats
#$ -M jennifer.semple@izb.unibe.ch
#$ -m ae
#$ -o /home/pmeister/log_files/
#$ -e /home/pmeister/error_logs/
#$ -t 1-2

#load R and libraries
module add R/3.2.2;
R -e 'library("genomation")'
R -e 'library("rtracklayer")'
R -e 'library("dplyr")'
R -e 'library("optparse")'

#RDIR=/home/pmeister/R-3.2.1/bin
	
DATADIR=/data1/projects/p025/Nano_C/Full_runs/FRG_Files

numColours=70
mkdir -p $DATADIR/bedFiles_${numColours}colours
mkdir -p $DATADIR/plots
mkdir -p $DATADIR/stats
mkdir -p $DATADIR/hops

#to process a batch of files using the ls command. You must set the -t option 
#above to go from 1 to the total number of files. These are then processed in parallel
files2process=(`ls $DATADIR/*LAST*.tsv`)
let i=$SGE_TASK_ID-1

Rscript bedFromSplitTSV.R -s -c $numColours ${files2process[$i]}
# -s means calculate fragment statistics
# -c $numColours gives the number of colours to use in bedfiles (between 12 and 70)



