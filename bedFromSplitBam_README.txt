#Basic workflow to create .bed files from .sam alignments:

# convert sam to bam 
samtools view -Sb  Meister_RCA1_N2_2_split_no_min_gm.sam >  Meister_RCA1_N2_2_split_no_min_gm.bam

# run python script bedFromSplitBam.py with --nsort option in order to sort bam by name (output is _nsorted.bam)
# and process data to output a .bed file
python bedFromSplitBam.py --nsort Meister_RCA1_N2_2_split_no_min_gm.bam

### the .bed file is good for viewing in text as you can see all the fragments from a single mother read consecutively 
### in the file. But to view in IGV tools you must sort again according to genomic position:

#for viewing in IGV sort the .bed file with IGVTools (can also be done from IGV in "Tools" menu)
#output is a _sorted.bed file
~/IGVTools/igvtools sort Meister_RCA1_N2_2_split_no_min_gm.bed Meister_RCA1_N2_2_split_no_min_gm_sorted.bed

#then index the file with IGVTools (also done automatically when you load into IGV)
#output is a .idx file
~/IGVTools/igvtools index Meister_RCA1_N2_2_split_no_min_gm_sorted.bed


