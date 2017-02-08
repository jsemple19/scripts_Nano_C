# 2017-20-01
# plotStatsXvsA.R

library(dplyr)


args <- commandArgs(TRUE)

for (fileName in args) {

  #setwd("~/Documents/MeisterLab/otherPeopleProjects/Nano_C/Frg_Files")
  #fileName="./hops/TSV_NoAmp2-Pass_SpltAln_LAST-Trained_Flattened.csv"
  dataDir<-dirname(fileName)
  #setwd(dataDir)
  fileBase<-strsplit(basename(fileName),".",fixed=TRUE)[[1]][1]
  hops<-read.csv(fileName,stringsAsFactors=FALSE)
  
  hops$distanceClass<-factor(hops$distanceClass,levels=c("<100b","<1k","<10k","<100k","<1m","<10m","Far cis","Trans"))
  hops$interChr<-factor(hops$interChr)
  #names(levels(hops$interChr))<-c("intraChr","interChr")
  
  ############## funcitons ##############
  doBarplot<-function(myData,...){
    x<-barplot(myData,ylab="counts",xlab="sequential hop distance categories",
               ylim=c(0,max(myData)*1.1),cex.main=0.9,cex.sub=0.8,...)
    text(x,myData*1.05,myData)
  }
  
  doFreqBarplot<-function(myData,...){
    myData<-round(myData/sum(myData),2)
    x<-barplot(myData, ylab="frequency",xlab="sequential hop distance categories",
               ylim=c(0,max(myData)*1.1),cex.main=0.9,cex.sub=0.8,...)
    text(x,myData*1.05,myData)
  }
  
  getFragsPerRead<-function(myData) {
    myData$hopToFrom<-gsub("([0-9])_FNID","\\1__FNID",myData$hopToFrom,  perl=TRUE)
    ID1<-sapply(strsplit(myData$hopToFrom,"__",fixed=TRUE),'[[',1)
    ID2<-sapply(strsplit(myData$hopToFrom,"__",fixed=TRUE),'[[',2)
    ldata<-cbind(ID1,myData[,c("lchr","SNID")])
    rdata<-cbind(ID2,myData[,c("rchr","SNID")])
    names(ldata)<-c("FNID","chr","SNID")
    names(rdata)<-c("FNID","chr","SNID")
    dataByFrag<-unique(rbind(ldata,rdata))
    fragNum<-dataByFrag %>% count(SNID,sort=TRUE)
    fragPR<-data.frame(table(fragNum$n),stringsAsFactors=FALSE)
    names(fragPR)<-c("fragPerRead","freq")
    fragPR$fragPerRead<-as.numeric(fragPR$fragPerRead)
    return(fragPR)
  }
  
  ###############3 end of functions #####################33
  
  pdf(file=paste0(dataDir,"/../plots/XvA_",fileBase,".pdf"),title=fileBase,width=8,height=11,paper="a4")
  
  ################ do this hop centered
  
  #first take only intra chromosomal (purest comparison for distances)
  layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE))
  X<-which(hops$lchr=="chrX" & hops$rchr=="chrX")
  A<-which(hops$lchr==hops$rchr)
  xTitle="ChrX: only X-X hops"
  aTitle="ChrA: only intra A-A hops "
  hist(hops$distance[X],main=xTitle,sub=fileBase,col="purple",breaks=25,cex.sub=0.9,cex.main=0.9)
  hist(hops$distance[A],main=aTitle,sub=fileBase, col="gray",breaks=25,cex.sub=0.9,cex.main=0.9)
  
  doBarplot(table(hops$distanceClass[X]),main=xTitle,sub=fileBase,col="purple")
  doBarplot(table(hops$distanceClass[A]),main=aTitle,sub=fileBase,col="gray")
  
  doFreqBarplot(table(hops$distanceClass[X]),main=xTitle,sub=fileBase,col="purple")
  doFreqBarplot(table(hops$distanceClass[A]),main=aTitle,sub=fileBase,col="gray")
  
  ####  distances between any involving X  vs all autosomes #### not well normalized #########3
  X<-which(hops$lchr=="chrX" | hops$rchr=="chrX")
  A<-which(hops$lchr!="chrX" | hops$rchr!="chrX")
  xTitle="ChrX: X-X, X-A, A-X hops"
  aTitle="ChrA: intra A-A & inter A-A hops"
  #layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE))
  hist(hops$distance[X],main=xTitle,sub=fileBase,col="purple",breaks=25,cex.sub=0.9,cex.main=0.9)
  hist(hops$distance[A],main=aTitle,sub=fileBase, col="gray",breaks=25,cex.sub=0.9,cex.main=0.9)
  
  doBarplot(table(hops$distanceClass[X]),main=xTitle,sub=fileBase,col="purple")
  doBarplot(table(hops$distanceClass[A]),main=aTitle,sub=fileBase,col="gray")
  
  doFreqBarplot(table(hops$distanceClass[X]),main=xTitle,sub=fileBase,col="purple")
  doFreqBarplot(table(hops$distanceClass[A]),main=aTitle,sub=fileBase,col="gray")
  
  #### distances between chrl involving X  vs any autosome
  X<-which(hops$lchr=="chrX")
  A<-which(hops$lchr!="chrX")
  xTitle="ChrX: X-X, X-A hops"
  aTitle="ChrA: intra and inter A-A, and A-X hops"
  #layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE))
  hist(hops$distance[X],main=xTitle,sub=fileBase,col="purple",breaks=25,cex.sub=0.9,cex.main=0.9)
  hist(hops$distance[A],main=aTitle,sub=fileBase, col="gray",breaks=25,cex.sub=0.9,cex.main=0.9)
  
  doBarplot(table(hops$distanceClass[X]),main=xTitle,sub=fileBase,col="purple")
  doBarplot(table(hops$distanceClass[A]),main=aTitle,sub=fileBase,col="gray")
  
  doFreqBarplot(table(hops$distanceClass[X]),main=xTitle,sub=fileBase,col="purple")
  doFreqBarplot(table(hops$distanceClass[A]),main=aTitle,sub=fileBase,col="gray")
  
  
  ########### do it centered on mother reads
  ## take any mother read that involves an chrX vs all those that don't
  X<-which(hops$lchr=="chrX" | hops$rchr=="chrX")
  
  Xhops<-hops[hops$SNID %in% hops$SNID[X],]
  Ahops<-hops[!(hops$SNID %in% hops$SNID[X]),]
  xTitle="ChrX: Any mother read containing a chrX frag"
  aTitle="ChrA: Any mother read not containing a chrX frag"
  
  layout(matrix(c(1,1,1,2,3,4,5,5,5,6,6,6),4,3,byrow=TRUE))
  
  # fragments per read
  xFragPR<-getFragsPerRead(Xhops)
  x<-barplot(xFragPR$freq,names.arg=xFragPR$fragPerRead, xlab="Fragments per read",
             ylab="Count",main=xTitle,sub=fileBase,ylim=c(0,1.1*max(xFragPR$freq)),
             cex.main=0.8,cex.sub=0.8,col="purple")
  text(x,xFragPR$freq+0.05*max(xFragPR$freq),labels=xFragPR$freq)
  
  #plot intra vs inter chromosomal frequency
  intChr<-table(Xhops$interChr)
  names(intChr)<-c("intraChr","interChr")
  x<-barplot(intChr,ylab="counts",main=xTitle,sub=fileBase,
             ylim=c(0,max(intChr)*1.1),cex.main=0.8,cex.sub=0.8,col="purple")
  text(x,intChr*1.05,labels=intChr)
  
  intChr_freq=round(intChr/sum(intChr),2)
  x<-barplot(intChr_freq,ylab="frequency",main=xTitle,sub=fileBase,
             ylim=c(0,max(intChr_freq)*1.1),cex.main=0.8,cex.sub=0.8,col="purple")
  text(x,intChr_freq*1.05,labels=intChr_freq)
  
  #plot histogram of hop distances 
  Xhops$distance<-as.numeric(Xhops$distance)   
  options(scipen=1)
  hist(Xhops$distance,breaks=25,xlab="sequential hop distance",main=xTitle,
       sub=fileBase,cex.main=0.8,cex.sub=0.8,col="purple")
  
  # barplot of hop distances by category
  doBarplot(table(Xhops$distanceClass),main=xTitle,sub=fileBase,col="purple")
  doFreqBarplot(table(Xhops$distanceClass),main=xTitle,sub=fileBase,col="purple")
  
  ## same for mother frags with autosomes
  # fragments per read
  aFragPR<-getFragsPerRead(Ahops)
  x<-barplot(aFragPR$freq,names.arg=aFragPR$fragPerRead, xlab="Fragments per read",
             ylab="Count",main=aTitle,sub=fileBase,ylim=c(0,1.1*max(aFragPR$freq)),
             cex.main=0.8,cex.sub=0.8,col="gray")
  text(x,aFragPR$freq+0.05*max(aFragPR$freq),labels=aFragPR$freq)
  
  #plot intra vs inter chromosomal frequency
  intChr<-table(Ahops$interChr)
  names(intChr)<-c("intraChr","interChr")
  x<-barplot(intChr,ylab="counts",main=aTitle,sub=fileBase,
             ylim=c(0,max(intChr)*1.1),cex.main=0.8,cex.sub=0.8)
  text(x,intChr*1.05,labels=intChr)
  
  intChr_freq=round(intChr/sum(intChr),2)
  x<-barplot(intChr_freq,ylab="frequency",main=aTitle,sub=fileBase,
             ylim=c(0,max(intChr_freq)*1.1),cex.main=0.8,cex.sub=0.8)
  text(x,intChr_freq*1.05,labels=intChr_freq)
  
  #plot histogram of hop distances 
  Ahops$distance<-as.numeric(Ahops$distance)   
  options(scipen=1)
  hist(Ahops$distance,breaks=25,xlab="sequential hop distance",main=aTitle,
       sub=fileBase,cex.main=0.8,cex.sub=0.8)
  
  # barplot of hop distances by category
  doBarplot(table(Ahops$distanceClass),main=aTitle,sub=fileBase,col="gray")
  doFreqBarplot(table(Ahops$distanceClass),main=aTitle,sub=fileBase,col="gray")
  
  
  
  ###### look at mother reads which are only intrachromosomal
  # take only intraChr mother reads
  A<-which(hops$lchr!="chrX" | hops$rchr!="chrX")
  X<-which(hops$lchr!=hops$rchr | hops$lchr=="chrX" | hops$rchr=="chrX")
  Xhops<-hops[!(hops$SNID %in% hops$SNID[A]),]
  Ahops<-hops[!(hops$SNID %in% hops$SNID[X]),]
  xTitle="ChrX: only intrachromosomal mother reads"
  aTitle="ChrA: only intrachromosomal mother reads"
  
  layout(matrix(c(1,1,1,2,3,4,5,5,5,6,6,6),4,3,byrow=TRUE))
  
  # fragments per read
  xFragPR<-getFragsPerRead(Xhops)
  x<-barplot(xFragPR$freq,names.arg=xFragPR$fragPerRead, xlab="Fragments per read",
             ylab="Count",main=xTitle,sub=fileBase,ylim=c(0,1.1*max(xFragPR$freq)),
             cex.main=0.8,cex.sub=0.8,col="purple")
  text(x,xFragPR$freq+0.05*max(xFragPR$freq),labels=xFragPR$freq)
  
  #plot intra vs inter chromosomal frequency
  intChr<-table(Xhops$interChr)
  names(intChr)<-c("intraChr","interChr")
  x<-barplot(intChr,ylab="counts",main=xTitle,sub=fileBase,
             ylim=c(0,max(intChr)*1.1),cex.main=0.8,cex.sub=0.8,col="purple")
  text(x,intChr*1.05,labels=intChr)
  
  intChr_freq=round(intChr/sum(intChr),2)
  x<-barplot(intChr_freq,ylab="frequency",main=xTitle,sub=fileBase,
             ylim=c(0,max(intChr_freq)*1.1),cex.main=0.8,cex.sub=0.8,col="purple")
  text(x,intChr_freq*1.05,labels=intChr_freq)
  
  #plot histogram of hop distances 
  Xhops$distance<-as.numeric(Xhops$distance)   
  options(scipen=1)
  hist(Xhops$distance,breaks=25,xlab="sequential hop distance",main=xTitle,
       sub=fileBase,cex.main=0.8,cex.sub=0.8,col="purple")
  
  # barplot of hop distances by category
  doBarplot(table(Xhops$distanceClass),main=xTitle,sub=fileBase,col="purple")
  doFreqBarplot(table(Xhops$distanceClass),main=xTitle,sub=fileBase,col="purple")
  
  ## same for mother frags with only autosomes
  # fragments per read
  aFragPR<-getFragsPerRead(Ahops)
  x<-barplot(aFragPR$freq,names.arg=aFragPR$fragPerRead, xlab="Fragments per read",
             ylab="Count",main=aTitle,sub=fileBase,ylim=c(0,1.1*max(aFragPR$freq)),
             cex.main=0.8,cex.sub=0.8,col="gray")
  text(x,aFragPR$freq+0.05*max(aFragPR$freq),labels=aFragPR$freq)
  
  #plot intra vs inter chromosomal frequency
  intChr<-table(Ahops$interChr)
  names(intChr)<-c("intraChr","interChr")
  x<-barplot(intChr,ylab="counts",main=aTitle,sub=fileBase,
             ylim=c(0,max(intChr)*1.1),cex.main=0.8,cex.sub=0.8)
  text(x,intChr*1.05,labels=intChr)
  
  intChr_freq=round(intChr/sum(intChr),2)
  x<-barplot(intChr_freq,ylab="frequency",main=aTitle,sub=fileBase,
             ylim=c(0,max(intChr_freq)*1.1),cex.main=0.8,cex.sub=0.8)
  text(x,intChr_freq*1.05,labels=intChr_freq)
  
  #plot histogram of hop distances 
  Ahops$distance<-as.numeric(Ahops$distance)   
  options(scipen=1)
  hist(Ahops$distance,breaks=25,xlab="sequential hop distance",main=aTitle,
       sub=fileBase,cex.main=0.8,cex.sub=0.8)
  
  # barplot of hop distances by category
  doBarplot(table(Ahops$distanceClass),main=aTitle,sub=fileBase,col="gray")
  doFreqBarplot(table(Ahops$distanceClass),main=aTitle,sub=fileBase,col="gray")
  
  
  dev.off()
}