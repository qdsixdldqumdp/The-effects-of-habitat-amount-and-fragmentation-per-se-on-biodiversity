setwd("/home/zoushuai/hasize/0")
.libPaths("/home/zoushuai/R/x86_64-pc-linux-gnu-library/4.0")
library(NLMR)
library(raster)
library(virtualspecies)
library(reshape2)
library(ggplot2)
library(ncdf4)
library(dplyr)
options(scipen = 10)


#envir,整数
inisimu<-function(nspn,nmapx,nmapy,envir){
  spn<<-nspn;mapx<<-nmapx;mapy<<-nmapy;mspp<-list()
  for (n in 1:spn) {
    mspp1 <- nlm_gaussianfield(mapx,mapy,resolution =1, autocorr_range =envir,mag_var = 1, nug = 0.1, mean = 0.5, rescale = TRUE)
    tranmspp<-convertToPA(mspp1,plot = F,alpha =-0.1,beta=0.5)
    mspp2<-as.matrix(tranmspp$probability.of.occurrence)
    mspp[[n]]<-mspp2
  }
  msp<<-mspp
  promsp<<-mspp
  
  #生成虚拟物种实况底图
  maplen<<-nmapx*nmapy
  themapp<<-matrix(nrow=mapy,ncol=mapx)
  #生成虚拟物种分布，物种存在由适合度决定。
  
  
  
  
  
  themap<<-themapp
  spdismap<<-themapp
}#初始化 给出变量物种数spn、地图x轴长mapx、地图y轴长mapy，给出物种适应性地图msp/promsp，分布图themap/spdismap
#------------------------



#disper :near radius global
simu<-function(times,disper,ppool){
  tim<-0
  qaq<-matrix(nrow = times,ncol = spn+1)
  colnames(qaq)<-0:spn
  while (tim<times) {
    tim<-tim+1
    
    for(xi in 1:mapx){
      for (yj in 1:mapy) {
        if(!is.na(themap[yj,xi])){
          idk<-themap[yj,xi]
          if(idk>0){
            suit<-msp[[idk]][yj,xi]
            if(rbinom(1,1,suit)!=1){
              if(rbinom(1,1,ppool)==1){
                poolsp<-sample(1:spn,1)
                themap[yj,xi]<-poolsp
              }else{
                themap[yj,xi]<-NA
                if(disper=="near"){
                  fi<-xi+c(-1,1);fii<-fi[fi>0&fi<=mapx]#列
                  fj<-yj+c(-1,1);fjj<-fj[fj>0&fj<=mapy]#行，并剔除超界
                  fec_sp<-NULL
                  for (feci in fii) {
                    if(!is.na(themap[yj,feci])){
                      fec_sp<-c(fec_sp,themap[yj,feci])
                    }}
                  for (fecj in fjj) {
                    if(!is.na(themap[fecj,xi])){
                      fec_sp<-c(fec_sp,themap[fecj,xi])
                    }}
                  sppool<-as.data.frame(table(fec_sp))
                  sppoolvec<-NULL
                  for(sppoolind in 1:length(sppool[,1])){
                    spind<-as.numeric(as.character(sppool[sppoolind,1]))
                    sppool[sppoolind,3]<-msp[[spind]][yj,xi]*sppool[sppoolind,2]
                    sppoolvec<-c(sppoolvec,spind)
                  }
                  if(length(sppoolvec)==0){
                    print("error")
                  }else if(length(sppoolvec)==1){
                    themap[yj,xi]<-sppoolvec
                  }else{  
                    sppoolprob<-as.numeric(sppool[,3])
                    themap[yj,xi]<-sample(sppoolvec,size=1,prob = sppoolprob)
                  } 
                }else if(disper=="radius"){
                  fi<-xi+c(-5:5);fii<-fi[fi>0&fi<=mapx]#列
                  fj<-yj+c(-5:5);fjj<-fj[fj>0&fj<=mapy]#行，并剔除超界
                  fec_sp<-NULL
                  for (feci in fii) {
                    for (fecj in fjj) {
                      if(!is.na(themap[fecj,feci])){
                        fec_sp<-c(fec_sp,themap[fecj,feci])
                      }}}
                  sppool<-as.data.frame(table(fec_sp))
                  sppoolvec<-NULL
                  for(sppoolind in 1:length(sppool[,1])){
                    spind<-as.numeric(as.character(sppool[sppoolind,1]))
                    sppool[sppoolind,3]<-msp[[spind]][yj,xi]*sppool[sppoolind,2]
                    sppoolvec<-c(sppoolvec,spind)
                  }
                  if(length(sppoolvec)==0){
                    print("error")
                  }else if(length(sppoolvec)==1){
                    themap[yj,xi]<-sppoolvec
                  }else{  
                    sppoolprob<-as.numeric(sppool[,3])
                    themap[yj,xi]<-sample(sppoolvec,size=1,prob = sppoolprob)
                  } 
                  
                }else if(disper=="global"){
                  sppool<-as.data.frame(table(themap))
                  sppoolvec<-NULL
                  for(sppoolind in 1:length(sppool[,1])){
                    spind<-as.numeric(as.character(sppool[sppoolind,1]))
                    sppool[sppoolind,3]<-msp[[spind]][yj,xi]*sppool[sppoolind,2]
                    sppoolvec<-c(sppoolvec,spind)
                  }
                  if(length(sppoolvec)==0){
                    print("error")
                  }else if(length(sppoolvec)==1){
                    themap[yj,xi]<-sppoolvec
                  }else{  
                    sppoolprob<-as.numeric(sppool[,3])
                    themap[yj,xi]<-sample(sppoolvec,size=1,prob = sppoolprob)
                  } 
                }else{print("erro")}
                
              }
            }
          }
        }
      }
    }
    
    mat<-themap
    for (mi in 0:spn) {qaq[tim,mi+1]<-length(which(mat==mi))}
  

  hax1<-7
  hax2<-18
  hay1<-7
  hay2<-18
  mat<-themap[hay1:hay2,hax1:hax2]
  for (mi in 0:spn) {qaq[tim,mi+1]<-length(which(mat==mi))}
  }
  pho<<-qaq;spdismap<<-themap
}#调用themap，msp。给出记录物种数目变化的矩阵pho,物种分布图spdismap

#enviro#空间自相关值 0-1
#disper#扩散距离 0-mapx/y
#pool#从物种库引入的概率





path<-
  
for(replic in 1:50){
  for(enviro in c(1,3,6,9,12)){
    for(pool in c(0.0001,0.001,0.01,0.1,1)){
      for (disper in c("near","radius","global")) {
        inisimu(80,32,32,enviro)
        
        simu(10000,dis=disper,ppool = pool)
        write.csv(pho,paste("rep",replic,"piex",path,disper,"env",enviro,"pool",pool,".csv"),row.names=F)
      }
    }
  }
}
