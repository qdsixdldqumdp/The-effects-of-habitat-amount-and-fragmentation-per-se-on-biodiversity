abundata<-read.csv(file.choose())
metadata<-read.csv(file.choose())
head(metadata)
slossmeta<-metadata[,c(2,13,15)]
getwd()
setwd("/Users/jingwufeng/Downloads/SLOSS")
library(ggplot2)
library(openxlsx)



studyind<-1
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
 slosstest1[i,7]<-slossmeta[pan1[pan2],3]
}
omitcon<-slosstest1[,7]
omitind<-which("continuous"==omitcon)
slosstest1[,7][omitind]<-max(as.numeric(omitcon[-omitind]))*3



                    
sloss1<-slosstest1[order(as.numeric(slosstest1$size)),]
sitepool<-names(table(slosstest1$site_id))
cunchu<-matrix(ncol=4,nrow=length(sitepool),0)
rownames(cunchu)<-paste(sitepool,"simulation")
colnames(cunchu)<-c("sl","sp1","ls","sp2")
n<-1;nn<-1
sppool<-NULL
while(n<=length(sloss1[,1])){ 
  sitename<-sloss1[n,3]
  siteind<-which(sitename==sloss1[,3])

  for(sitei in 1:length(siteind)){
    spname<-sloss1[siteind[sitei],5]
    spsimulation<-length(which(spname==sppool))
    if(spsimulation==0){sppool<-c(sppool,spname)}
  }
  cunchu[nn,2]<-as.numeric(length(sppool))
  if(nn<=1){cunchu[nn,1]<-as.numeric(sloss1[n,7])
  }else{
    nm<-nn-1
    cunchu[nn,1]<-as.numeric(sloss1[n,7])+cunchu[nm,1]
    }
  nn<-nn+1
  n<-siteind[length(siteind)]+1
}

  
  
sloss1<-slosstest1[order(as.numeric(slosstest1$size),decreasing = TRUE),]

  n<-1;nn<-1
  sppool<-NULL
  while(n<=length(sloss1[,1])){ 
    sitename<-sloss1[n,3]
    siteind<-which(sitename==sloss1[,3])
    
    for(sitei in 1:length(siteind)){
      spname<-sloss1[siteind[sitei],5]
      spsimulation<-length(which(spname==sppool))
      if(spsimulation==0){sppool<-c(sppool,spname)}
    }
    cunchu[nn,4]<-as.numeric(length(sppool))
    if(nn<=1){cunchu[nn,3]<-as.numeric(sloss1[n,7])
    }else{
      nm<-nn-1
      cunchu[nn,3]<-as.numeric(sloss1[n,7])+cunchu[nm,3]
    }
    nn<-nn+1
    n<-siteind[length(siteind)]+1
  }

  p<-ggplot(data = as.data.frame(cunchu))+
    geom_point(aes(x=sl,y=sp1,color="red"))+
    geom_point(aes(x=ls,y=sp2,color="blue"))+
    geom_line(aes(x=sl,y=sp1,color="red"))+
    geom_line(aes(x=ls,y=sp2,color="blue"))+
    labs(x="area simulation",y="sp simulation",title = myfilename)+
    scale_color_discrete(name="sp-area curve",breaks=c("red","blue"),
                         labels=c("Small to Large","Large to Small"))+
    theme_classic()
  p
    
ggsave(paste(myfilename,".pdf"),p,device = "pdf")
write.csv(cunchu,paste(myfilename,".csv"))
}
  
