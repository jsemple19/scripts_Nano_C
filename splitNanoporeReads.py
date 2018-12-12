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
.fastq extension. Option to give output file names does not work, but you can specify
an output directory. Default is ../NanoC_out/splitFastq

python ./splitNanoporeReads.py input1.fastq input2.fastq input3.fastq --outputFileDir outputDir

The script outputs fastq files split into fragments at GATC sites as independent reads
It also outputs a seqData.txt with a list of all the fragments being created along
with the read data and run data.
"""
import re
import argparse
import ntpath
import os

#oldfastq="/home/jenny/Documents/Meister_RCA1_N2_2cp.fastq"
#newfastq="/home/jenny/Documents/Meister_RCA1_N2_split.fastq"
#minFragLen=20

def splitNanoporeFastQ(oldfastq,newfastq,minFragLen):
    # open files for reading and writing
    with open(oldfastq, 'r') as oldfq, open(newfastq, 'w') as newfq, open(seqDataFile,'a') as seqData:
        while True:
            # start reading lines and search for one staring with @ (ID line of read)
            uniqueID=oldfq.readline()
            if not uniqueID: 
                break
            if not uniqueID.startswith("@"):
                continue
            seq=oldfq.readline() # read in sequence line
            line3=oldfq.readline() # read in "+" line 3
            qual=oldfq.readline() # read in quality line
            totalReadLength=len(seq) # get total read length
            if not qual:
                break
            # find all GATCs in sequence
            starts=[m.start() for m in re.finditer("GATC",seq)]
            begin=0
            # extract the read and run ID from name
            readID=re.sub("^@","",uniqueID.split()[0])
            runID=re.sub("runid=","",uniqueID.split()[1])
            #cycle through the GATC positions
            for i in range(len(starts)):
                end=starts[i]
                # create unique id for the fragment with ordinal number and start:end
                fragID=readID+"_frag"+"{0:0>3}".format(str(i+1))+"_"+str(begin+1)+":"+str(end)
                # write each fragment that is greater thab minSize as a separate read (4 lines)
                if((end-begin)>minFragLen):
                    newfq.write("@"+fragID+"\n")
                    newfq.write(seq[begin:end]+"\n")
                    newfq.write(line3)
                    newfq.write(qual[begin:end]+"\n")
                # write read and frag ids and lengths to seqData.txt file
                seqData.write(readID+"\t"+str(totalReadLength)+"\t"+fragID+"\t"+str(end-begin)+"\t"+runID+"\n")
                begin=end
            # do the same for last fragment in read
            fragID=readID+"_frag"+"{0:0>3}".format(str(len(starts)+1))+"_"+str(begin+1)+":"+str(len(seq))
            if((len(seq)-begin)>minFragLen):
                end=len(seq)
                newfq.write("@"+fragID+"\n")
                newfq.write(seq[begin:-1]+"\n")
                newfq.write(line3)
                newfq.write(qual[begin:-1]+"\n")
            seqData.write(readID+"\t"+str(totalReadLength)+"\t"+fragID+"\t"+str(end-begin)+"\t"+runID+"\n")
                

#splitNanoporeFastQ(oldfastq,newfastq,minFragLen)              

if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Split NanoPore reads at GATC")
    parser.add_argument('inputFile',nargs='*',help="Name of fastq format input file(s)")
    parser.add_argument('--outputFile',help="Name of output file (optional)")
    parser.add_argument('--outputDir',help="Name of output directory (for multi-file processing)",default="../NanoC_out")
    parser.add_argument('--minFragLen',nargs='?',default=20,type=int,help="Discard fragments shorter than this length (default=20)")

    args=parser.parse_args()
    minFragLen=args.minFragLen #smallest fragment size that would be of interest for mapping
    outputDir=args.outputDir
    seqDataFile=outputDir+"/seqData.txt"
    seqData=open(seqDataFile,'w')
    seqData.write("readID\treadLength\tfragmentID\tfragmentLength\tfastqID\n")
    seqData.close()
    for oldfastq in args.inputFile:
    #oldfastq=args.inputFile
        if ((args.outputFile==None) | (len(args.inputFile)>1)):
            if not os.path.exists(outputDir+"/splitFastq"):
                os.makedirs(outputDir+"/splitFastq")
            newfastq=outputDir+"/splitFastq/"+re.sub("\.fastq$","_split.fastq",ntpath.basename(oldfastq))
        else:
            outputDir="."
            newfastq=args.outputFile
        splitNanoporeFastQ(oldfastq,newfastq,minFragLen)

    
    
    
