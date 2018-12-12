### Values that need changing
setwd("~/Documents/MeisterLab/otherPeopleProjects/Moushumi/13102018_hic2/scripts_Nano_C")
expName="13102018_hic2"
path="../NanoC_out/"
importantStats<-list()

### Initiation variables, packages, etc
library(lattice)
require(stats)

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)

#Chrnames<- c("CHROMOSOME_I","CHROMOSOME_II","CHROMOSOME_III","CHROMOSOME_IV","CHROMOSOME_V","CHROMOSOME_X")
Chrnames<-c("I","II","III","IV","V","X") # used to get rid of mtDNA
chrOrder<-c("I","II","III","IV","V","X","MtDNA") # used to reorder sequence name levels in objects


### read in data
fragAlin<-read.delim(paste0(path,"aln/", expName,"_q10.tab"),header=T,stringsAsFactors=F)
fragID<-read.table(paste0(path,"seqData.txt"),header=T,stringsAsFactors=F,sep="\t")
fragID$readID<-gsub("@","",fragID$readID)
names(FragID)<-c("readFileIDs","fastqFile","headerData")
IDs<-strsplit(FragID$readFileIDs,";")

uReadID<-sapply(IDs,"[[",3)
FragID$ReadId<-as.numeric(gsub("uRN:","",uReadID))


FragAlin1<-merge(FragAlin,FragID,by="ReadId")
write.csv(FragAlin1,paste0(path, expName,"_withReadID.csv"),row.names=FALSE)

importantStats["numAlignedFrag"]<-dim(FragAlin1)[1]
importantStats["numReadsWithAlignedFrag"]<-length(unique(FragAlin1$ReadId))

##############################################################################
### Basic Stats from HiC Nanopore Datasets - After pymcHiC pipeline
#############################################################################
# Adriana's plots
#################  Read file

FragAlin<-read.csv(paste0(path, expName,"_withReadID.csv"),header=T,stringsAsFactors=F)
FragAlin <- FragAlin[FragAlin$AlnChr %in% Chrnames,]
#################

### Addition of Fragment Length to table and change of Class of ReadId
FragLength <- FragAlin$FragmentEndBp - FragAlin$FragmentStartBp
FragAlin$FragmentLen <- FragLength
#str(FragAlin) # To get info
FragAlin$ReadId <- as.factor(FragAlin$ReadId)

importantStats["medianReadLength"]<-median(FragAlin$ReadLen)
importantStats["medianFragLength"]<-median(FragAlin$FragmentLen)
importantStats["minReadLength"]<-min(FragAlin$ReadLen)
importantStats["minFragLength"]<-min(FragAlin$FragmentLen)
importantStats["maxReadLength"]<-max(FragAlin$ReadLen)
importantStats["maxFragLength"]<-max(FragAlin$FragmentLen)

### Basic Stats start HERE
### Number of Fragments per Chr - NOTE:Some Fragments are duplicated
FragmChr <- table(FragAlin$AlnChr)
FragmChr <- FragmChr[Chrnames]
pdf(paste0(path, expName,"_FragAperChr.pdf"))  # Change name accordingly
par(las=2)
par(mar=c(10, 5, 4, 4))
par(mfrow=c(1,1))
barplot(FragmChr, main="#Fragments per Chromosome", col="lightgreen")
dev.off()


### Size distribution of Fragments per Chr - NOTE:Some Fragments are duplicated
x<- split(FragAlin$FragmentLen,FragAlin$AlnChr)
x<- x[Chrnames]
pdf(paste0(path, expName,"_FSizeperChr.pdf"))
par(las=2)
par(mar=c(10, 4, 4, 2))
par(mfrow=c(1,2))
boxplot(x,outline=FALSE, main="Fragment Size Distribution\nper Chr w/o Outliers", ylab="Size Distribution", col="lightgreen",boxwex = 0.25)
boxplot(x,main='Fragment Size Distribution\nper Chr', ylab="Size Distribution", col="lightgreen",boxwex = 0.25)
dev.off()

### Number of aligned reads(totally or partially)
y <- split(FragAlin$FragmentId,FragAlin$ReadId)
length(y)
### Fragments(non-duplicated) per read
d<-sapply(y,unique)
t<-sapply(d,length)
w <- sapply(t, function(h) if(h>=20){h<-20}else{h<-h})
w<-as.factor(w)
FragAlignxRead <- table(w)
write.table(FragAlignxRead, 'CircleSize.txt')
pdf(paste0(path, expName,"_CircleSize.pdf"))
par(las=2)
par(mar=c(4, 4, 4, 2))
barplot(FragAlignxRead, main="# of Reads containing\n# Fragments(Aligned)", col="lightgreen", xlab="# Fragment")
dev.off()

### Read size arranged by # of Fragments
w <- sapply(t, function(h) if(h>=20){h<-20}else{h<-h})
q <- split(FragAlin$ReadLen,FragAlin$ReadId)
q <-sapply(q,unique)
m <- data.frame(x=w,y=q)
m$x <- as.factor(m$x)
RLeng_CircleS <- split(m$y,m$x)
pdf(paste0(path, expName,"_ReadLeng_CircSize_1.pdf"))
par(las=2)
par(mar=c(4, 4, 4, 2))
boxplot(RLeng_CircleS,outline=FALSE, main="Read Size\nper #Fragment w/o Outliers", col="lightgreen",boxwex = 0.25)
dev.off()

pdf(paste0(path, expName,"_ReadLeng_CircSize_2.pdf"))
par(las=2)
par(mar=c(4, 4, 4, 2))
boxplot(RLeng_CircleS,main='Read Size\nper #Fragment', ylab="Size Distribution", col="lightgreen",boxwex = 0.25)
dev.off()


#read with fragments
par(mfrow=c(1,1))





##############################################################################
### summarise by read, the number of fragments mapped and the number of chormosomes
#############################################################################
# Jenny's old plots

options(scipen=10^9)

i=match("barcode",names(FragAlin))
if(is.na(i)){
  FragAlin$barcode<-0
}
readData<-FragAlin %>%
  group_by(barcode,ReadId) %>%
  summarise(
    numFrags=length(unique(FragmentId)),
    numChr=length(unique(AlnChr)))

#summarise counts by barcode for plotting
dt<-readData %>%
  group_by(barcode,numFrags) %>%
  summarise(
    count=n())

p1<-ggplot(dt,aes(x=as.factor(numFrags),y=count)) +
  geom_bar(stat="identity") +
  facet_wrap(~barcode) + ggtitle("Number of fragments per read") +
  scale_x_discrete(name="Fragments per read")


#summarise frequency by barcode for plotting
dt<-readData %>%
  group_by(barcode,numFrags) %>%
  summarise(count=n()) %>%
  mutate(frequency=count/sum(count))

p2<-ggplot(dt,aes(x=as.factor(numFrags),y=frequency,fill=as.factor(barcode))) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c('red', 'blue')) + labs(fill="barcode") +
  ggtitle("Frequency of N fragments per read") +
  scale_x_discrete(name="Fragments per read (N)")

# plot a closeup of 1-10 fragments
p3<-ggplot(dt,aes(x=as.factor(numFrags),y=frequency,fill=as.factor(barcode))) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c('red', 'blue')) + labs(fill="barcode") +
  ggtitle("Frequency of N fragments per read (up to 10 frags only)") +
  scale_x_discrete(name="Fragments per read (N)",limits=c(1:10))

pdf(paste0(path, expName,"_fragmentsPerRead.pdf"),paper="a4",height=11,width=8)
grid.arrange(p1,p2,p3,nrow=3)
dev.off()

#filter for reads with 2 or more fragments
multiFrag<-readData[readData$numFrags>1,]
importantStats["numMultiFragReads"]<-length(unique(multiFrag$ReadId))

## plot counts of number of chormosomes per read.

#summarise counts by barcode for plotting
dt<-multiFrag %>%
  group_by(barcode,numChr) %>%
  summarise(
    count=n())

p4<-ggplot(dt,aes(x=as.factor(numChr),y=count)) +
  geom_bar(stat="identity") +
  facet_wrap(~barcode) + ggtitle("Chromosomes per read of multi-fragment reads") +
  scale_x_discrete(name="Chromosomes per read")


#summarise frequency by barcode for plotting
dt<-multiFrag %>%
  group_by(barcode,numChr) %>%
  summarise(count=n()) %>%
  mutate(frequency=count/sum(count))

p5<-ggplot(dt,aes(x=as.factor(numChr),y=frequency,fill=as.factor(barcode))) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c('red', 'blue')) + labs(fill="barcode") +
  ggtitle("Frequency of multi-fragment reads with N chromosomes per read") +
  scale_x_discrete(name="Chromosomes per read (N)")

pdf(paste0(path, expName,"_chrsPerRead.pdf"),paper="a4",height=11,width=8)
grid.arrange(p4,p5,nrow=3)
dev.off()

# filter for reads that hit multiple chromosomes and save as .csv
multiChr<-multiFrag[multiFrag$numChr>1,]
idx<-FragAlin$ReadId %in% multiChr$ReadId
write.csv(FragAlin[idx,],paste0(path, expName,"_readsWithMultipleChrHits.csv"),row.names=F)




##############################################################################
### make table of chromosomal contacts (hops)
#############################################################################

#####################
##### functions
####################


minDistance<-function(v1,v2,v3,v4) {
  # function finds the smallest distance between all fragment ends
  #it takes four vectors representing AlnStart, AlnEnd, lead(AlnStart), lead(AlnEnd)
  # output is a vector with minimum absolute distance
  allDist<-cbind(v1-v3,v2-v4,v2-v3,v1-v4)
  minDist<-apply(abs(allDist),1,min)
  return(minDist)
}

pairOrientation<-function(df) {
  # function finds the relative orientation of two fragments
  # input is a dataframe with five columns: AlnStart,AlnEnd, lead(AlnStart),lead(AlnEnd) AlnStrand
  # output is one of "divergent","convergent","tandem","inverted".
  df1<-df
  # reverse orientation of start end according to strand
  idx<-which(df[,5]=="False")
  df1[idx,1]<-df[idx,2]
  df1[idx,2]<-df[idx,1]
  idx<-which(df[,6]=="False")
  df1[idx,3]<-df[idx,4]
  df1[idx,4]<-df[idx,3]
  df<-df1
  allDist<-cbind(as.numeric(df[,1])-as.numeric(df[,3]),
                 as.numeric(df[,2])-as.numeric(df[,4]),
                 as.numeric(df[,2])-as.numeric(df[,3]),
                 as.numeric(df[,1])-as.numeric(df[,4]))
  # find rows with NA
  naRows<-(rowSums(is.na(allDist))>0)
  minType<-apply(abs(allDist[!naRows,]),1,which.min)
  ori<-c("divergent","convergent","tandem","inverted")[minType]
  allOri<-rep(NA,length(naRows))
  allOri[!naRows]<-ori
  return(allOri)
}


isOverlapping<-function(df) {
  # function finds fragment pairs from contacts that overlap
  # output is TRUE (overlapping) or FALSE (not overlapping)
  GR1<-GRanges(seqnames=df[,5],ranges=IRanges(start=as.numeric(df[,1]),end=as.numeric(df[,2])), strand="*")
  GR2<-GRanges(seqnames=df[-dim(df)[1],6],ranges=IRanges(start=as.numeric(df[-dim(df)[1],3]),
                                              end=as.numeric(df[-dim(df)[1],4])), strand="*")
  ol<-findOverlaps(GR1,GR2)
  idx<-which(queryHits(ol)==subjectHits(ol))
  overlapping<-rep(FALSE, length(GR1))
  overlapping[queryHits(ol)[idx]]<-TRUE
  return(overlapping)
}


getOverlapSize<-function(df) {
  # function finds fragment pairs from contacts that overlap
  # output is TRUE (overlapping) or FALSE (not overlapping)
  GR1<-GRanges(seqnames=df[,5],ranges=IRanges(start=as.numeric(df[,1]),end=as.numeric(df[,2])), strand="*")
  GR2<-GRanges(seqnames=df[-dim(df)[1],6],ranges=IRanges(start=as.numeric(df[-dim(df)[1],3]),
                                                         end=as.numeric(df[-dim(df)[1],4])), strand="*")
  ol<-findOverlaps(GR1,GR2)
  idx<-which(queryHits(ol)==subjectHits(ol))
  sizeOverlap<-rep(0, length(GR1))
  sizeOverlap[queryHits(ol)[idx]]<-width(pintersect(GR1[queryHits(ol)[idx]],GR2[subjectHits(ol)[idx]]))
  return(sizeOverlap)
}


fixOverlapOri<-function(df) {
  # this function fixes the orientation calculated by pairOrientation which can fail if fragments
  # overlap.
  # output is vector with orientations
  names(df)<-c("overlap","sameStrand","AlnStrand","AlnStart1","AlnStart2","hopOri","ReadId")
  hopOri<-df$hopOri
  df<-df[,-6]
  i<-which(df$overlap==T & df$sameStrand==T & df$AlnStrand=="True")
  df1<-df[i,]
  df1$hopOri<-with(df1,
    ifelse(AlnStart1<AlnStart2,"tandem","inverted"))
  hopOri[i]<-df1$hopOri
  ##
  i<-which(df$overlap==T & df$sameStrand==T & df$AlnStrand=="False")
  df1<-df[i,]
  df1$hopOri<-with(df1,
                   ifelse(AlnStart1>AlnStart2,"tandem","inverted"))
  hopOri[i]<-df1$hopOri
  ##
  i<-which(df$overlap==T & df$sameStrand==F & df$AlnStrand=="True")
  df1<-df[i,]
  df1$hopOri<-with(df1,
                   ifelse(AlnStart1<AlnStart2,"convergent","divergent"))
  hopOri[i]<-df1$hopOri
  ##
  i<-which(df$overlap==T & df$sameStrand==F & df$AlnStrand=="False")
  df1<-df[i,]
  df1$hopOri<-with(df1,
                   ifelse(AlnStart1>AlnStart2,"convergent","divergent"))
  hopOri[i]<-df1$hopOri
  ##
  return(hopOri)
}

###############
##############
###############

### contact table

# contact table has several additional fields as follows:
# AlnChr2 - chromosome to which the second fragment aligns to
# hopType - within same chromosome (Intra) or between chromosomes (Inter)
# hopSize - minimal absolute distance between fragment ends
# sameStrand - are both fragments on the both strand?
# hopOri - orientation of hop:
#       tandem:     -frag1-->  -frag2-->  or  <--frag2- <--frag1-
#       inverted:   -frag2--> -frag1-->   or  <--frag1- <--frag2-
#       convergent: -frag1--> <--frag2-   or  -frag2--> <--frag1-
#       divergent:  <--frag2- -frag1-->   or  <--frag1- -frag2-->
# hopOverlap - do the aligned fragments overlap?
# firstChr - chr,strand,start and end of alignment of first fragment in pair
# secondChr - chr, strand, start and end of alignment of second fragment in pair
FragAlin<-read.csv(paste0(path, expName,"_withReadID.csv"),header=T,stringsAsFactors=F)

# remove fragments that align more than one place
FragAlin$uniqueID<-unite(FragAlin,"uniqueID",c("ReadId","FragmentId"))$uniqueID
FragAlin<-FragAlin[!duplicated(FragAlin$uniqueID),]

# change default display of tibble to show all columns
options(tibble.width = Inf)

# hops<-FragAlin %>%
#   distinct(ReadId, FragmentId, .keep_all = TRUE) %>% # remove duplicated fragments
#   arrange(ReadId,FragmentId) %>% # sort by read id then fragment id
#   mutate(hopOri=pairOrientation(cbind(AlnStart,AlnEnd,lead(AlnStart),lead(AlnEnd),AlnStrand,lead(AlnStrand)))) %>%
#   mutate(overlap=isOverlapping(cbind(AlnStart,AlnEnd,lead(AlnStart),lead(AlnEnd),AlnChr,lead(AlnChr)))) %>%
#   group_by(ReadId) %>%
#   mutate(AlnChr2=lead(AlnChr)) %>%
#   mutate(hopType=ifelse(FragmentId!=lead(FragmentId) & AlnChr==lead(AlnChr),"Intra","Inter")) %>%
#   mutate(hopSize=minDistance(AlnStart,AlnEnd,lead(AlnStart),lead(AlnEnd))) %>%
#   mutate(sameStrand=ifelse(AlnStrand==lead(AlnStrand),TRUE,FALSE)) %>%
#   mutate(firstChr=paste0(AlnChr,ifelse(AlnStrand=="True","+","-"),":",AlnStart,"-",AlnEnd)) %>%
#   mutate(secondChr=paste0(lead(AlnChr),ifelse(lead(AlnStrand)=="True","+","-"),":",
#                           lead(AlnStart),"-",lead(AlnEnd)))


# make a table of combinatorial combinations of fragments
multiFrag<-FragAlin %>%
  arrange(ReadId,FragmentId) %>% # sort by read id then fragment id
  dplyr::group_by(ReadId) %>%
  dplyr::summarise(fragNum=n()) %>%
  dplyr::filter(fragNum>1)

FragAlin<-FragAlin[FragAlin$ReadId %in% multiFrag$ReadId,]

cols2combine<-FragAlin[,c("ReadId","FragmentId")] %>%
  dplyr::group_by(ReadId) %>%
  tidyr::expand(nesting(ReadId,FragmentId),nesting(ReadId,FragmentId)) %>%
  filter(ReadId==ReadId1 & FragmentId<FragmentId1) # remove self interactions and duplicate two way interactions

cols2combine<-cols2combine %>%
  mutate(consecutive=ifelse((FragmentId1-FragmentId)==1,T,F)) # check if fragments in pair are immediately consecutive to each other

# combine ReadId and FragmentId to get unique ID for both first and second fragment
cols2combine$uniqueID1<-tidyr::unite(cols2combine,"uniqueID1",c("ReadId","FragmentId"))$uniqueID1
cols2combine$uniqueID2<-tidyr::unite(cols2combine,"uniqueID2",c("ReadId1","FragmentId1"))$uniqueID2

# add second fragment data for table
cols2combine<-merge(FragAlin,cols2combine[,c("consecutive","uniqueID1","uniqueID2")],by.x="uniqueID",by.y="uniqueID1")
hopsComb<-cols2combine
idx<-match(hopsComb$uniqueID2,FragAlin$uniqueID)
hopsComb$AlnChr2<-FragAlin[idx,"AlnChr"]
hopsComb$AlnStart2<-FragAlin[idx,"AlnStart"]
hopsComb$AlnEnd2<-FragAlin[idx,"AlnEnd"]
hopsComb$AlnStrand2<-FragAlin[idx,"AlnStrand"]


hopsComb<-hopsComb %>%
  mutate(hopOri=pairOrientation(cbind(AlnStart,AlnEnd,AlnStart2,AlnEnd2,AlnStrand,AlnStrand2))) %>%
  mutate(overlap=isOverlapping(cbind(AlnStart,AlnEnd,AlnStart2,AlnEnd2,AlnChr,AlnChr2))) %>%
  mutate(sizeOverlap=getOverlapSize(cbind(AlnStart,AlnEnd,AlnStart2,AlnEnd2,AlnChr,AlnChr2))) %>%
  mutate(hopType=ifelse(AlnChr==AlnChr2,"Intra","Inter")) %>%
  mutate(hopSize=minDistance(AlnStart,AlnEnd,AlnStart2,AlnEnd2)) %>%
  mutate(sameStrand=ifelse(AlnStrand==AlnStrand2,TRUE,FALSE)) %>%
  mutate(firstChr=paste0(AlnChr,ifelse(AlnStrand=="True","+","-"),":",AlnStart,"-",AlnEnd)) %>%
  mutate(secondChr=paste0(AlnChr2,ifelse(AlnStrand2=="True","+","-"),":",AlnStart2,"-",AlnEnd2))


# when hops are interchromosomal, the stand, size and orientation of the hop are meaningless
hopsComb$hopSize[hopsComb$hopType=="Inter"]<-NA
hopsComb$sameStrand[hopsComb$hopType=="Inter"]<-NA
hopsComb$hopOri[hopsComb$hopType=="Inter"]<-NA


# correct hop orientation when fragments overlap
df<-hopsComb %>%
  dplyr::select(overlap,sameStrand,AlnStrand,AlnStart,AlnStart2,hopOri,ReadId)
hopsComb$hopOri<-fixOverlapOri(df)

hopsComb<-hopsComb %>%
  arrange(ReadId,FragmentId)

# save table
write.csv(hopsComb,paste0(path,expName,"_contacts.csv"),row.names=F)




##############################################################################
### plot basic contact stats
#############################################################################

contacts<-read.csv(paste0(path,expName,"_contacts.csv"),stringsAsFactors=F)

# overlappingFrags<-contacts %>%
#   filter(overlap==TRUE) %>%
#   select(-c("fastqFile","headerData"))
# #46878 overlapping from 346339 contacts
#
# write.csv(overlappingFrags,paste0(path,expName,"_overlappingFrags.csv"),row.names=F)


p1<-ggplot(contacts,aes(x=hopType,fill=hopType)) + geom_bar() + scale_fill_brewer(palette="Dark2") +
  ggtitle("Inter- vs intra-chromosomal contacts")

p2<-contacts %>%
  filter(!is.na(hopOri)) %>%
  ggplot(aes(x=hopOri,fill=hopOri))  + geom_bar() + scale_fill_brewer(palette="Dark2") +
  ggtitle("Orientation of fragments")

p3<-contacts %>%
  filter(!is.na(hopOri)) %>%
  ggplot(aes(x=overlap,fill=overlap))  + geom_bar() + facet_wrap(~hopOri) +
  scale_fill_brewer(palette="Dark2") +
  ggtitle("Overlap of fragments")

p4<- contacts %>%
  filter(!is.na(hopOri)) %>%
  ggplot(aes(x=sameStrand,fill=sameStrand))  + geom_bar() + facet_wrap(~hopOri) +
  scale_fill_brewer(palette="Dark2") +
  ggtitle("Strandedness of fragments")

p5<- contacts %>%
  filter(!is.na(hopOri)) %>%
  ggplot(aes(x=sameStrand,y=hopSize,fill=sameStrand))  + geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = quantile(contacts$hopSize,na.rm=T, c(0.1, 0.9))) +
  scale_fill_brewer(palette="Dark2") +
  ggtitle("Distance between fragments")

p6<- contacts %>%
  filter(!is.na(hopOri)) %>%
  ggplot(aes(x=hopOri,y=hopSize,fill=hopOri))  + geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = quantile(contacts$hopSize,na.rm=T, c(0.1, 0.9))) +
  scale_fill_brewer(palette="Dark2") +
  ggtitle("Distance between fragments")

p7 <- contacts %>%
  filter(!is.na(hopOri)) %>%
  ggplot(aes(x=hopSize))  + geom_histogram(bins=60) +
  ggtitle("Distance between fragments")

minHop=500
p8 <- contacts %>%
  filter(!is.na(hopOri)) %>%
  filter(hopSize>minHop) %>%
  ggplot(aes(x=hopSize))  + geom_histogram(bins=60) +
  ggtitle(paste0("Distance between fragments (>",minHop,"bp)"))

p9 <- contacts %>%
  filter(!is.na(hopOri)) %>%
  ggplot(aes(x=hopSize,fill=sameStrand))  + geom_histogram(bins=60) + facet_wrap(~sameStrand) +
  scale_fill_brewer(palette="Dark2") +
  ggtitle(paste0("Distance between fragments by strandedness"))

p10 <- contacts %>%
  filter(!is.na(hopOri)) %>%
  ggplot(aes(x=hopSize,fill=hopOri))  + geom_histogram(bins=60) + facet_wrap(~hopOri) +
  scale_fill_brewer(palette="Dark2") +
  ggtitle(paste0("Distance between fragments by orientation"))

ml<-marrangeGrob(grobs=list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10),ncol=2,nrow=1)
ggsave(paste0(path, expName,"_basicContactStats.pdf"),ml,device="pdf",width=29,height=20,
       units="cm")

######### look at data at the per read level

contacts<-contacts %>%
  group_by(ReadId) %>%
  mutate(numContacts=n())

# compare number of inter vs intra chromosomal contacts by number of contacts in a read
dt<-contacts %>%
  group_by(ReadId,numContacts,hopType) %>%
  summarise(count=n())

p1<-ggplot(dt,aes(x=as.factor(numContacts),y=count,fill=as.factor(hopType))) +
      geom_boxplot()  +
      ggtitle("Proportion of interchromosomal contacts vs number of contacts") +
  scale_x_discrete(name="Contacts per read")

# compare number of inter vs intra chromosomal contacts by read length
dt<-contacts %>%
  group_by(ReadId,ReadLen,hopType) %>%
  summarise(count=n())

p2<-ggplot(dt,aes(x=ReadLen,y=count,col=as.factor(hopType))) +
  geom_point()  + geom_smooth() +
  ggtitle("Proportion of interchromosomal contacts vs read length") +
  scale_x_discrete(name="Read length")

# compare number of inter vs intra chromosomal contacts per read
dt<-contacts %>%
  group_by(ReadId,ReadLen,hopType) %>%
  summarise(count=n()) %>%
  spread(key=hopType,value=count,fill=0)

p3<-ggplot(dt,aes(x=Intra,y=Inter)) +
  geom_jitter(col="darkblue",alpha=0.3) +
  ggtitle("Inter vs Intra contacts per read") +
  scale_x_continuous(name="Number of Intra-chromosomal contacts per read") +
  scale_y_continuous(name="Number of Inter-chromosomal contacts per read")

ml<-marrangeGrob(grobs=list(p1,p2,p3),ncol=1,nrow=1)
ggsave(paste0(path, expName,"_contactsByRead.pdf"),ml,device="pdf",width=20,height=20,
       units="cm")

# compare number of intra and inter chromosomal contacts by chr type


# first lets plot interchr contacts by individual chromosomes, and use the first chr
# to determine chr type (This is not very good becuase for interchromsomal contacts, second
# chromosome could be chrX)

contacts$AlnChr<-factor(contacts$AlnChr,levels=chrOrder)
contacts$AlnChr2<-factor(contacts$AlnChr2,levels=chrOrder)

dt<-contacts %>%
  group_by(AlnChr,hopType) %>%
  summarise(count=n()) %>%
  mutate(frequency=count/sum(count))

xlevel<-as.numeric(dt[dt$AlnChr=="X" & dt$hopType=="Intra","frequency"])
p4<-ggplot(dt,aes(x=AlnChr,y=frequency,fill=hopType)) + geom_bar(stat="identity") +
  geom_hline(yintercept=xlevel,linetype="dashed",
             color = "grey", size=0.5) +
  ggtitle("Intra vs inter contacts by chr (!!!!!!only chr1 considered)")

## if just looking at Intra chromsomal contacts, looking at first fragment chr is ok

dt<-contacts %>%
  filter(hopType=="Intra") %>%
  group_by(AlnChr) %>%
  summarise(count=n())

p5<-ggplot(dt,aes(x=AlnChr,y=count)) + geom_bar(stat="identity") +
  geom_hline(yintercept=xlevel,linetype="dashed",
             color = "grey", size=0.5) +
  ggtitle("Intra chromosomal contacts by chr")


dt$normContacts<-dt$count/seqlengths(Celegans)
xlevel<-as.numeric(dt[dt$AlnChr=="X","normContacts"])
p6<-ggplot(dt,aes(x=AlnChr,y=normContacts)) + geom_bar(stat="identity") +
  geom_hline(yintercept=xlevel,linetype="dashed",
             color = "grey", size=0.5) +
  ggtitle("Intra chromosomal contacts by chr, normalised by chr length")


dt<-dt %>%
  filter(dt$AlnChr!="MtDNA")
xlevel<-as.numeric(dt[dt$AlnChr=="X","normContacts"])
p7<-ggplot(dt,aes(x=AlnChr,y=normContacts)) + geom_bar(stat="identity") +
  geom_hline(yintercept=xlevel,linetype="dashed",
             color = "grey", size=0.5) +
  ggtitle("Intra chromosomal contacts by chr, normalised by chr length")



# for inter chromosomal contact we will define ANY contact that involves an X chr, as an X chr contact
anyX<-grepl("X",contacts$firstChr) | grepl("X",contacts$secondChr)
contacts$chrType<-ifelse(anyX,"chrX","autosomes")
# extract name of second Chr in pair
#contacts$AlnChr2<-gsub("[[:punct:]]?:[[:digit:]]*-[[:digit:]]*$","",contacts$secondChr)
#contacts$AlnChr2<-factor(contacts$AlnChr2,levels=chrOrder)

dt<-contacts %>%
  filter(hopType=="Inter") %>%
  group_by(AlnChr,AlnChr2) %>%
  summarise(count=n())

p8<-ggplot(dt,aes(x=AlnChr2,y=count)) + geom_bar(stat="identity") +
  facet_wrap(~AlnChr) +
  ggtitle("Inter chromosomal contacts (panel=frag1,bar=frag2)")

# normalise to chr length
seqnames(Celegans)<-chrOrder
idx<-match(dt$AlnChr2,seqnames(Celegans))
dt$normContacts<-dt$count/seqlengths(Celegans)[idx]

p9<-ggplot(dt,aes(x=AlnChr2,y=normContacts)) + geom_bar(stat="identity") +
  facet_wrap(~AlnChr) +
  ggtitle("Normalised inter chromosomal contacts (panel=frag1,bar=frag2)")

# remove mtDNA
dt<- dt %>%
  filter(AlnChr!="MtDNA" & AlnChr2!="MtDNA")

p10<-ggplot(dt,aes(x=AlnChr2,y=count)) + geom_bar(stat="identity") +
  facet_wrap(~AlnChr) +
  ggtitle("Inter chromosomal contacts (panel=frag1,bar=frag2)")


p11<-ggplot(dt,aes(x=AlnChr2,y=normContacts)) + geom_bar(stat="identity") +
  facet_wrap(~AlnChr) +
  ggtitle("Normalised inter chromosomal contacts (panel=frag1,bar=frag2)")


idx<-match(dt$AlnChr,seqnames(Celegans))
dt$normContacts2x<-dt$count/seqlengths(Celegans)[idx]

p12<-ggplot(dt,aes(x=AlnChr2,y=normContacts2x)) + geom_bar(stat="identity") +
  facet_wrap(~AlnChr) +
  ggtitle("Normalised (x2) inter chromosomal contacts (panel=frag1,bar=frag2)")

ml<-marrangeGrob(grobs=list(p4,p5,p6,p7,p8,p9,p10,p11,p12),ncol=1,nrow=1)
ggsave(paste0(path, expName,"_contactsByChr.pdf"),ml,device="pdf",width=20,height=20,
       units="cm")

