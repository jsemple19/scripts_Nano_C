
"""
2018-12-11
script to create dataframe from mapped split reads in a bam file.
The script extracts basic read data such as 
readId: fragment name 
chr: chromosome to which it aligns 
alnStart,alnEnd: start, end (calculated by removing softclipping from ref start end
strand: strand to which it aligns
mapQ: mapping quality score 
refStart, refEnd: reference start, end using pysam function 
queryStart, quaryEnd: query start, end using pysam function
fragLength: length of fragment
Qscore from bam file and saves it as a .tab file.
"""

import pysam
import re
import argparse
import ntpath



def dfFromSplitBam(bamFileName):
    # extract the output directory from the input file
    outDir=ntpath.dirname(bamFileName)
    # extract the data set name from the input file.
    datasetName=ntpath.basename(bamFileName).split(".bam")[0]
    datasetName=re.sub("_nsorted","",datasetName)
    # read in the bam file
    myBam=pysam.AlignmentFile(bamFileName,"rb")
    # open the output file and write the header to it
    outTab=open(outDir+"/"+datasetName+".tab",'w')
    header=["readId","chr","alnStart","alnEnd","strand","mapQ","refStart","refEnd","queryStart","queryEnd","fragLength"]
    outTab.write("\t".join(header)+"\n")
    for seg in myBam:
        # only look at mapped sequences
        if seg.is_unmapped==False:
            # extract strand info
            if seg.is_reverse==True:
                strand="-"
            else:
                strand="+"
            # extract cigar info to calculate match start and end
            cigar=seg.cigartuples
            # check for soft clipping on the left hand side
            if cigar[0][0]==4:
                LF=cigar[0][1]
            else:
                LF=0
            # check for soft clipping on the right hand side
            if cigar[-1][0]==4:
                RF=cigar[-1][1]
            else:
                RF=0
            # get inferred fragment length for comparison
            fragLen=seg.infer_query_length()
            # get reference start and end to calculate match start and end
            refStart=seg.reference_start
            refEnd=seg.reference_end
            # get fragment ID
            ID=seg.query_name
            # get mapping quality
            mapQ=seg.mapping_quality
            # get query start and end for comparison (aligned portion of sequence - anything that is not soft clipped)
            queryStart=seg.query_alignment_start
            queryEnd=seg.query_alignment_end
            # get chr name to which it aligns
            chrm=myBam.get_reference_name(seg.reference_id)
            # write data for this alignment to output file
            tabLine=[ID,chrm,str(refStart-LF),str(refEnd+RF),strand,str(mapQ),str(refStart),str(refEnd),str(queryStart),str(queryEnd),str(fragLen)]
            outTab.write("\t".join(tabLine)+"\n")  


if __name__=='__main__':
    parser=argparse.ArgumentParser(description="create bedfile for split NanoPore reads")
    parser.add_argument('inputFile',help="Name of .bam input file")
    args=parser.parse_args()
    bamFileName=args.inputFile
    dfFromSplitBam(bamFileName)
