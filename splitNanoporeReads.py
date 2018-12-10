#! /usr/bin/python
"""
splitNanoporeReads.py
2017-01-03
Script to read in a fastq file and split the read (and quality score) at GATC
sequence and create a new fastq file with read name preserved
usage:

1) Process single file,output is stored in input_split.fastq (modification of
input file name). All fragments shorter than 20 nt (default) are discarded: 

python ./splitNanoporeReads.py input.fastq

2)Process single file and output is stored in output.fastq (name provided by 
user).All fragments shorter than 50 nt are discarded:

python ./splitNanoporeReads.py input.fastq --outputFile output.fastq --minFragLen 50

3) Process multiple files (output file names created by inserting _split before
.fastq extension. option to give output file names does not work)

python ./splitNanoporeReads.py input1.fastq input2.fastq input3.fastq
"""
import re
import argparse
#oldfastq="/home/jenny/Documents/Meister_RCA1_N2_2cp.fastq"
#newfastq="/home/jenny/Documents/Meister_RCA1_N2_split.fastq"
#minFragLen=20

def splitNanoporeFastQ(oldfastq,newfastq,minFragLen):
    with open(oldfastq, 'r') as oldfq, open(newfastq, 'w') as newfq open(seqData,'r') as seqData:
        seqData.write("readID\treadLength\tfragmentID\tfragmentLength\n")
        while True:
            uniqueID=oldfq.readline()
            if not uniqueID: 
                break
            if not uniqueID.startswith("@"):
                continue
            seq=oldfq.readline()
            line3=oldfq.readline()
            qual=oldfq.readline()
            totalReadLength=len(seq)
            if not qual:
                break
            starts=[m.start() for m in re.finditer("GATC",seq)]
            begin=0
            fragIDbase=re.sub("Basecall_Alignment_template","",uniqueID.split()[0])
            for i in range(len(starts)):
                end=starts[i]
                fragID=fragIDbase+"_frag"+"{0:0>3}".format(str(i+1))+"_"+str(begin+1)+":"+str(end)
                if((end-begin)>minFragLen):
                    newfq.write(fragID+"\n")
                    newfq.write(seq[begin:end]+"\n")
                    newfq.write(line3)
                    newfq.write(qual[begin:end]+"\n")
                seqData.write(fragIDbase+"\t"+str(totalReadLength)+"\t"+fragID+"\t"+str(end-begin+"\n")
                begin=end
            fragID=fragIDbase+"frag"+"{0:0>3}".format(str(len(starts)+1))+"_"+str(begin+1)+":"+str(len(seq))
            if((len(seq)-begin)>minFragLen):
                newfq.write(fragID+"\n")
                newfq.write(seq[begin:-1]+"\n")
                newfq.write(line3)
                newfq.write(qual[begin:-1]+"\n")

                

#splitNanoporeFastQ(oldfastq,newfastq,minFragLen)              

if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Split NanoPore reads at GATC")
    parser.add_argument('inputFile',nargs='*',help="Name of fastq format input file(s)")
    parser.add_argument('--outputFile',help="Name of output file (optional)")
    parser.add_argument('--minFragLen',nargs='?',default=20,type=int,help="Discard fragments shorter than this length (default=20)")

    args=parser.parse_args()
    minFragLen=args.minFragLen #smallest fragment size that would be of interest for mapping
    for oldfastq in args.inputFile:
    #oldfastq=args.inputFile
        if ((args.outputFile==None) | (len(args.inputFile)>1)):
            newfastq=re.sub("\.fastq$","_split.fastq",oldfastq)
        else:
            newfastq=args.outputFile
        seqData="./seqData.txt"
        splitNanoporeFastQ(oldfastq,newfastq,minFragLen)

    
    
    
