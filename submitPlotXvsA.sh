#!/bin/bash
#$ -l h_cpu=72:00:00
#$ -l h_vmem=8G
#$ -cwd
#$ -N XvsA_stats
#$ -M jennifer.semple@izb.unibe.ch
#$ -m ae
#$ -o /home/pmeister/log_files/
#$ -e /home/pmeister/error_logs/


#load R and libraries
module add R/3.2.2;
R -e 'library("dplyr")'

#RDIR=/home/pmeister/R-3.2.1/bin
	
DATADIR=/data1/projects/p025/Nano_C/Full_runs/FRG_Files

mkdir -p $DATADIR/plots

#to process a batch of files using the ls command. You must set the -t option 
#above to go from 1 to the total number of files. These are then processed in parallel
files2process=(`ls $DATADIR/hops/*Pass*LAST*.csv`)

for filename in "${files2process[@]}"
do
	Rscript plotStatsXvsA.R "$filename"
done

