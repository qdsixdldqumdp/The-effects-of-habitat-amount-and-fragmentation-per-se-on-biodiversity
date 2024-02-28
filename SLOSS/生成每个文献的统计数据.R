abundata<-read.csv(file.choose())
metadata<-read.csv(file.choose())
head(metadata)
slossmeta<-metadata[which("balanced"==metadata$sampling_effort_unit),c(2,13,15)]
getwd()
setwd("/Users/jingwufeng/Downloads/habitat fragment/SLOSS/equal-size&landscape/data")
library(ggplot2)
library(openxlsx)



studyind<-1
author1<-names(table(slossmeta$refshort))
metastudy<-matrix(nrow = length(author1),ncol = 4)
metastudy[,1]<-author1
colnames(metastudy)<-c("author","Intercept","slope","p_slope")


while(studyind<length(abundata$id)){
studyvec<-which(abundata[studyind,2]==abundata[,2])
myfilename<-abundata[studyind,2]
studyind<-1+studyvec[length(studyvec)]
slosstest1<-abundata[studyvec,]
colnames(slosstest1)[7]<-"size"
slosstest1[,7]<-"size"
for (i in 1:length(studyvec)) {
  pan1<-which(slosstest1[i,2]==slossmeta[,1])
  
 pan2<-which(slosstest1[i,3]==slossmeta[pan1,2])
 if(length(pan2)>1){
 if(slossmeta[pan1[pan2][1],3]==slossmeta[pan1[pan2][2],3]){
   pan2<-pan2[1]
 }}
 if(length(pan2)==1){slosstest1[i,7]<-slossmeta[pan1[pan2],3]}
}
omitcon<-slosstest1[,7]
omitind<-which("continuous"==omitcon)
slosstest1[,7][omitind]<-max(as.numeric(omitcon[-omitind]))*3

sloss1<-slosstest1[order(as.numeric(slosstest1$size)),]
cunchurow<-length(table(sloss1$site_id))
cunchu<-matrix(nrow = cunchurow,ncol = 4)
cunchu[,4]<-1:cunchurow
cunchu<-as.data.frame(cunchu)
colnames(cunchu)<-c("sp","size","plot name","x")
n<-1;nn<-1

while(n<=length(sloss1[,1])){ 
  sitename<-sloss1[n,3]
  cunchu[nn,3]<-sloss1[n,3]
  siteind<-which(sitename==sloss1[,3])
  cunchu[nn,1]<-length(table(sloss1[siteind,5]))
  cunchu[nn,2]<-sloss1[n,7]
  nn<-nn+1
  n<-siteind[length(siteind)]+1
}

s<-summary(lm(cunchu$sp~cunchu$x))
aa<-s$coefficients
indrow<-which(myfilename==metastudy[,1])
metastudy[indrow,2]<-round(aa[1,1],4)
metastudy[indrow,3]<-round(aa[2,1],4)
metastudy[indrow,4]<-round(aa[2,4],4)
write.csv(cunchu,paste(myfilename,".csv"))
}

setwd("/Users/jingwufeng/Downloads/habitat fragment/SLOSS/equal-size&landscape")

emoind<-which(NaN==metastudy[,4])
metastudy1<-metastudy[-emoind,]
write.csv(metastudy1,paste("metastudy",".csv"))

