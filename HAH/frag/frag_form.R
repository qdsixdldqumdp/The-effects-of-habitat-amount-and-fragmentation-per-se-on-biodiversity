setwd("/home/zoushuai/frag")
.libPaths("/home/zoushuai/R/x86_64-pc-linux-gnu-library/4.0")
library(NLMR)
library(raster)
library(virtualspecies)
library(reshape2)
library(stringr)
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
  themapp<-matrix(nrow=mapy,ncol=mapx)
  #生成虚拟物种分布，物种存在由适合度决定。
  themapp[1:maplen]<-sample(1:spn,maplen,replace = T)
  themap<<-themapp
  spdismap<<-themapp
}#初始化 给出变量物种数spn、地图x轴长mapx、地图y轴长mapy，给出物种适应性地图msp/promsp，分布图themap/spdismap
#------------------------#

map1patches<-list(x=c(11,22),y=c(11,22))
map2patches<-list(x=c(5,10,22,27),y=c(11,22))
map4patches<-list(x=c(5,10,22,27),y=c(5,10,22,27))
map8patches<-list(x=c(5,6,12,13,19,20,26,27),y=c(5,13,19,27))
map16patches<-list(x=c(4,6,11,13,18,20,25,27),y=c(4,6,11,13,18,20,25,27))
map36patches<-list(x=c(3,4,8,9,13,14,18,19,23,24,28,29),y=c(3,4,8,9,13,14,18,19,23,24,28,29))
map72patches<-list(x=c(0,1,4,5,8,9,12,13,16,17,20,21,24,25,27,28,30,31),y=c(2,2,6,6,10,10,14,14,18,18,22,22,26,26,30,30))
map144patches<-list(x=seq(0,35,by=3),y=seq(0,35,by=3))

fragmap<-function(nmap,fmapx=mapx,fmapy=mapy){
  fmap<-matrix(nrow = fmapx,ncol = fmapy)
  xyrange<-get(str_glue("map{nmap}patches"))
  findx<-length(xyrange$x)/2
  findy<-length(xyrange$y)/2
  for(indxx in 1:findx){
    for(indyy in 1:findy){
      indx<-indxx*2-1
      indy<-indyy*2-1
      x1<-xyrange$x[indx]+1
      x2<-xyrange$x[indx+1]+1
      y1<-xyrange$y[indy]+1
      y2<-xyrange$y[indy+1]+1
      fmap[y1:y2,x1:x2]<-themap[y1:y2,x1:x2]
    }
  }
  themap<<-fmap
  #plot(raster(themap))
}

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
  }
  pho<<-qaq;spdismap<<-themap
}#调用themap，msp。给出记录物种数目变化的矩阵pho,物种分布图spdismap

#enviro#空间自相关值 0-1
#disper#扩散距离 0-mapx/y
#pool#从物种库引入的概率

for(replic in 11:50){
  for(enviro in c(1,3,6,9,12)){
    for(pool in c(0.0001,0.001,0.01,0.1,1)){
      for (disper in c("near","radius","global")) {
        inisimu(80,36,36,enviro)
        
        fragmap(patchi)
        simu(10000,dis=disper,ppool = pool)
        write.csv(pho,paste("rep",replic,"map",mapx,disper,"env",enviro,"pool",pool,".csv"),row.names=F)
      }
    }
  }
}
