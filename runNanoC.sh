#! /bin/bash

expName="13102018_hic2"
fastqFiles=(../*_hic_*/fastq/pass/*)
genomeFile=~/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa
blackListFile=~/
outputDir=../NanoC_out
mkdir -p ${outputDir}
echo "splitting reads..."
python ./splitNanoporeReads.py ${fastqFiles[@]}

mergedFastq=${outputDir}/splitFastq/${expName}_merged.fastq

echo "merging all fastq files into one..."
cat ${outputDir}/splitFastq/*_split.fastq > ${mergedFastq}

rm  ${outputDir}/splitFastq/*_split.fastq

#if [[ ! -e "${outputDir}/aln" ]]
mkdir -p ${outputDir}/aln

echo "aligning to genome..."
bwa mem ${genomeFile} ${mergedFastq} > ${outputDir}/aln/${expName}.sam

# keep 0nly mapped reads that are primar alignments.
# Multiple mapping :One of these alignments is considered primary. 
# All the other alignments have the secondary alignment flag set in the SAM records that represent them.
# 256 for secondary alignment
# 4 for unmapped
# 2048 for supplementary alignment - do not remove this as they are chimeric alignments that could
# arise from missing a GATC due ot sequence errors. 
# (256 + 4  -> 260)
samtools view -F 260 ${outputDir}/aln/${expName}.sam -b -o ${outputDir}/aln/${expName}_q0.bam
samtools view -F 260 ${outputDir}/aln/${expName}.sam -q 10 -b -o ${outputDir}/aln/${expName}_q10.bam

# remove unnecessary sam file
#rm ${outputDir}/aln/${expName}.sam

#bedtools intersect -abam file.bam -b filter.bed -v > filtered.bam

# sort bam files by name
#samtools sort -n ${outputDir}/aln/${expName}_q0.bam -o ${outputDir}/aln/${expName}_q0_nsorted.bam
#samtools sort -n ${outputDir}/aln/${expName}_q10.bam -o ${outputDir}/aln/${expName}_q10_nsorted.bam

# extract data into a tab separated file 
python ./dfFromSplitBam.py ${outputDir}/aln/${expName}_q0.bam
python ./dfFromSplitBam.py ${outputDir}/aln/${expName}_q10.bam