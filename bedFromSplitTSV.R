#2017-01-19
#bedFromSplitTSV.R
#script to take TSV files supplied by Amin Allahyar of Nanopore split read mapping
#and convert it to a bed file with unique colours for each mother read

library(genomation)
library(rtracklayer)
library(dplyr)
library(optparse)

option_list <- list(
  make_option(c("-s", "--calcStats"), action="store_true", default=TRUE,
              help="calculate fragment statistics [default]"),
  make_option(c("-c", "--numColours"), type="integer", default=12,
              help="Number of colours to use in Bed file [default %default]",
              metavar="number")
)

#opt <- parse_args(OptionParser(option_list=option_list))
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = 1)
opt <- arguments$options
filename <- arguments$args
numColours<-opt$numColours
calcStats<-opt$calcStats

#filename= "/data1/projects/p025/Nano_C/Full_runs/FRG_Files/TSV_NoAmp1-Fail_SpltAln_LAST-Trained_Stitched.tsv"

#setwd("~/fromAminAllahyar")
#filename="./TSV_NoAmp1-Pass_SpltAln_LAST-Trained_Flattened.tsv"
#get commands from command line:
#args <- commandArgs(TRUE)


dataDir<-dirname(filename)
fileBase<-strsplit(basename(filename),".",fixed=TRUE)[[1]][1]
tsvData<-read.delim(filename,stringsAsFactors=FALSE)


############### some functions #####################
getSampleColours<-function(IDs2colour,numColours=12){
  # create different colour schemes. Not that 70 colours might overload the browser
  shades<-c(50,100,150,200)
  numShades<-length(shades)

  base6<-c("255,0,0","0,255,0","0,0,255","255,0,255","255,255,0","0,255,255")
  base12<-c(base6,"150,0,0","0,150,0","0,0,150","150,0,150","150,150,0","0,150,150")
  colExtend70<-c(paste(rep(shades,each=numShades*numShades),
                             rep(rep(shades,each=numShades),numShades),
                             rep(rep(shades,numShades),numShades),sep=","))
  shades<-c(50,125,200)
  numShades<-length(shades)
  colExtend27<-c(paste(rep(shades,each=numShades*numShades),
                       rep(rep(shades,each=numShades),numShades),
                       rep(rep(shades,numShades),numShades),sep=","))
  if (numColours<=12) {
    myColours<-base12
  } else if (numColours <= 27) {
    myColours<-c(base6,sample(colExtend27,numColours-6))
  } else if (numColours <= 70) {
    myColours<-c(base6,sample(colExtend70,numColours-6))
  } else {
    myColours<-c(base6,colExtend70)
  }
  return(sample(myColours,IDs2colour,replace=TRUE))
}


addLeading0s<-function(listNIDs) {
  decimalPlaces<-nchar(max(listNIDs))
  return(sprintf(paste0("%0",decimalPlaces,"d"),listNIDs))
}

################# end of functions ##########################

#randomly sample enough colours for all reads
sampleColour<-getSampleColours(max(tsvData$Sample_NID),numColours=numColours)

#add a sampleColour column to the .tsv file
tsvData$sampleColour<-sampleColour[tsvData$Sample_NID]
rm(sampleColour)

# add leading 0s onto sample ID for better sorting. Adjust the nubmer of 0s depending on the maximum sample number
tsvData$Sample_NID<-paste0("SNID_",addLeading0s(tsvData$Sample_NID),"_",tsvData$Read_Len)

# add leading 0s onto frag ID for better sorting. Adjust the nubmer of 0s depending on the maximum frag number
tsvData$Frg_NID<-  paste0("FNID_",addLeading0s(tsvData$Frg_NID),"_",tsvData$Seq_Be,"-",tsvData$Seq_En)

# change chromosome name to correct format
chr<-tsvData$Chr
chr<-factor(chr)
levels(chr)<-c(paste0("chr",c("I","II","III","IV","V","X","M")))
tsvData$Chr<-as.character(chr)

# change orientation to correct format
strand<-tsvData$Orientation
strand[strand==-1]<-"-"
strand[strand==0]<-"*"
strand[strand==1]<-"+"
tsvData$Orientation<-strand

#sort tsv file by sample and then fragment number
tsvData<-tsvData[order(tsvData$Sample_NID,tsvData$Frg_NID),]
tsvData$compositeID<-paste0(tsvData$Sample_NID,"__",tsvData$Frg_NID)

#find samples with more than 1 frag
nonUnique<-tsvData$Sample_NID[duplicated(tsvData$Sample_NID)]
tsvDataManyFrag<-tsvData[tsvData$Sample_NID %in% nonUnique,]

#create bed file for all data
bedFileName=paste0(dataDir,"/bedFiles_",numColours,"colours/",fileBase,".bed")
myBed<-tsvData[,c("Chr","Map_ExtBe","Map_ExtEn","compositeID","MQ","Orientation","Map_Be","Map_En","sampleColour")]
write.table(myBed,file=bedFileName,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

#create bed file for samples with multiple fragments only
bedFileNameMF=paste0(dataDir,"/bedFiles_",numColours,"colours/",fileBase,"_MF",".bed")
myBedMF<-tsvDataManyFrag[,c("Chr","Map_ExtBe","Map_ExtEn","compositeID","MQ","Orientation","Map_Be","Map_En","sampleColour")]
write.table(myBedMF,file=bedFileNameMF,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

################### calculate stats ####################
if (calcStats==TRUE) {
   rm(hops)
   curChr=''
   curSNID=''
   for (i in 1:dim(tsvDataManyFrag)[1]) {
      line<-tsvDataManyFrag[i,]
      if (curSNID!=line$Sample_NID | curSNID=='') {
         curSNID=line$Sample_NID
         curChr=line$Chr
         next
      } else {
         if (curChr==line$Chr) {
            interChr<-0
            lchr=curChr
            rchr=curChr
            if (line$Map_ExtBe>=tsvDataManyFrag[i-1,"Map_ExtEn"]) {
               distance<-line$Map_ExtBe-tsvDataManyFrag[i-1,"Map_ExtEn"]
               lpos=tsvDataManyFrag[i-1,"Map_ExtEn"]
               rpos=line$Map_ExtBe
            } else {
               distance<-tsvDataManyFrag[i-1,"Map_ExtBe"]-line$Map_ExtEn
               lpos=line$Map_ExtEn
               rpos=tsvDataManyFrag[i-1,"Map_ExtBe"]
            }
         } else {
            interChr<-1
            distance<-2.5e+07
            lchr=curChr
            rchr=line$Chr
            lpos=tsvDataManyFrag[i-1,"Map_ExtEn"]
            rpos=line$Map_ExtBe
            curChr<-line$Chr
         }
         hopFromTo<-paste0(tsvDataManyFrag[i-1,"Frg_NID"],"_",line$Frg_NID)
         if (!exists("hops")) {
            hops<-data.frame("hopToFrom"=hopFromTo,"lchr"=lchr,"lpos"=lpos,"rchr"=rchr,"rpos"=rpos,
                             "SNID"=curSNID,"interChr"=interChr,
                             "distance"=distance,stringsAsFactors=FALSE)
         } else {
            hops<-rbind(hops,c("hopToFrom"=hopFromTo,"lchr"=lchr,"lpos"=lpos,"rchr"=rchr,"rpos"=rpos,
                               "SNID"=curSNID,"interChr"=interChr,
                               "distance"=distance,stringsAsFactors=FALSE))
         }
      }
   }

   #calculate the number of fragments per read
   fragNum<-tsvData %>% count(Sample_NID,sort=TRUE)
   fragPR<-data.frame(table(fragNum$n),stringsAsFactors=FALSE)
   names(fragPR)<-c("fragPerRead","freq")
   fragPR$fragPerRead<-as.numeric(fragPR$fragPerRead)

   pdf(file=paste0(dataDir,"/plots/",fileBase,".pdf"),width=8,height=11,paper="a4")
   layout(matrix(c(1,1,1,2,3,4,5,5,5,6,6,6),4,3,byrow=TRUE))
   #plot barplot of the number of reads with x fragments
   x<-barplot(fragPR$freq,names.arg=fragPR$fragPerRead, xlab="Fragments per read",
              ylab="Count",main=fileBase,ylim=c(0,1.1*max(fragPR$freq)))
   text(x,fragPR$freq+0.05*max(fragPR$freq),labels=fragPR$freq)

   #plot intra vs inter chromosomal frequency
   intChr<-table(hops$interChr)
   names(intChr)<-c("intraChr","interChr")
   x<-barplot(intChr,ylab="counts",main=fileBase,ylim=c(0,max(intChr)*1.1))
   text(x,intChr*1.05,labels=intChr)

   intChr_freq=round(intChr/sum(intChr),2)
   x<-barplot(intChr_freq,ylab="frequency",main=fileBase,ylim=c(0,max(intChr_freq)*1.1))
   text(x,intChr_freq*1.05,labels=intChr_freq)

   #plot histogram of hop distances
   hops$distance<-as.numeric(hops$distance)
   options(scipen=1)
   hist(hops$distance,breaks=25,xlab="sequential hop distance",main=fileBase)

   # barplot of hop distances by category
   hops$distanceClass<-rep(NA,dim(hops)[1])
   hops$distanceClass[hops$distance>=1e+07]<-"Far cis"
   hops$distanceClass[hops$interChr==1]<-"Trans"
   hops$distanceClass[hops$distance<1e+07]<-"<10m"
   hops$distanceClass[hops$distance<1e+06]<-"<1m"
   hops$distanceClass[hops$distance<1e+05]<-"<100k"
   hops$distanceClass[hops$distance<1e+04]<-"<10k"
   hops$distanceClass[hops$distance<1e+03]<-"<1k"
   hops$distanceClass[hops$distance<1e+02]<-"<100b"
   hops$distanceClass<-factor(hops$distanceClass,levels=c("<100b","<1k","<10k","<100k","<1m","<10m","Far cis","Trans"))
   x<-barplot(table(hops$distanceClass), main=fileBase, ylab="counts",
              xlab="sequential hop distance categories",ylim=c(0,max(table(hops$distanceClass))*1.1))
   text(x,table(hops$distanceClass)*1.05,table(hops$distanceClass))

   distanceClass_freq<-round(table(hops$distanceClass)/sum(table(hops$distanceClass)),2)
   x<-barplot(distanceClass_freq, main=fileBase, ylab="frequency",
              xlab="sequencial hop distance categories",ylim=c(0,max(distanceClass_freq)*1.1))
   text(x,distanceClass_freq*1.05,distanceClass_freq)
   dev.off()

   # write summary of statistics and hop data to files
   statsSummary<-list("fragmentsPerRead"=fragPR,"intraVSinterChr"=intChr,"sequantialHopDistance"=table(hops$distanceClass))
   saveRDS(statsSummary,file=paste0(dataDir,"/stats/",fileBase,".Rds"))
   write.csv(hops,file=paste0(dataDir,"/hops/",fileBase,".csv"),row.names=FALSE,quote=FALSE)
}