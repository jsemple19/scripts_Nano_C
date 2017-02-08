
"""
2017-12-01
script to create bedfile from mapped reads in a bam file.
The .bam file must be sorted previously with the -n flag (sort by name). 
The script colours subfragments from the same mother fragment with the same colour
(randomly chosen from 9 different colours, so the same colour does not garauntee
that it is from the same mother read - double check fragment names).
The fragments are displayed as thick boxes representing the matched region and 
thin extensions which have the same length as the soft trimmed ends of the fragment
to give an idea of how far the fragment could reach in the best case scenario 
if its whole length was mapable.
To load the output file into IGV you must first sort it with IGV tools (tools 
menu in IGV), and create en index.
"""

import pysam
import random
import argparse

#bamFileName="Meister_RCA1_N2_2_split_no_min_bwasw_nsorted.bam"
#pysam.sort("-n",samFileName, "-o", datasetName+"_nsorted.bam")

def bedFromSplitBam(bamFileName,nsort=False):
    if nsort==True:
        datasetName=bamFileName.split(".bam")[0]
        pysam.sort("-n",bamFileName,"-o", datasetName+"_nsorted.bam")
        bamFileName=datasetName+"_nsorted.bam"
    datasetName=bamFileName.split("_nsorted.bam")[0]
    myBam=pysam.AlignmentFile(bamFileName,"rb")
    outBed=open(datasetName+".bed",'w')
    colours=["255,0,0","100,0,0","0,255,0","0,100,0","0,0,255","0,0,100","0,150,150","150,0,150","150,150,0"]
    IDcolour={}
    for seg in myBam:
        if seg.is_unmapped==False:
            if seg.is_reverse==True:
                strand="-"
            else:
                strand="+"
            cigar=seg.cigartuples
            if cigar[0][0]==4:
                LF=cigar[0][1]
            else:
                LF=0
            if cigar[-1][0]==4:
                RF=cigar[-1][1]
            else:
                RF=0
            #fragLen=seg.infer_query_length()
            refStart=seg.reference_start
            refEnd=seg.reference_end
            ID=seg.query_name
            motherID=ID.split("_")[0]
            if motherID not in IDcolour.keys():
                IDcolour[motherID]=random.sample(colours,1)
            #queryStart=seg.query_alignment_start
            #queryEnd=seg.query_alignment_end
            chrm=myBam.get_reference_name(seg.reference_id)
            bedLine=[chrm,str(refStart-LF-1),str(refEnd+RF-1),ID,"0",strand,str(refStart-1),str(refEnd-1),IDcolour[motherID][0]]
            outBed.write("\t".join(bedLine)+"\n")  


if __name__=='__main__':
    parser=argparse.ArgumentParser(description="create bedfile for split NanoPore reads")
    parser.add_argument('inputFile',help="Name of .bam input file")
    parser.add_argument('--nsort',action='store_true',help="sort the bam file by name")
    args=parser.parse_args()
    bamFileName=args.inputFile
    nsort=args.nsort
    bedFromSplitBam(bamFileName,nsort)
