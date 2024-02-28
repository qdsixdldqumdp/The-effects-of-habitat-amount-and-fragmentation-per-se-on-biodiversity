library(stringr)
studypool<-matrix(nrow=50*7*12*5*3,ncol = 8)
colnames(studypool)<-c("reptime","patch","enviro","poolposible","dispeisal","individual","sp","simpsen")
patchivec<-c(1, 2, 4, 8, 16, 36, 72,144)
indrow<-1
for(patchi in patchivec){
  setwd(str_glue("/home/zoushuai/frag/{patchi}"))
  for(replic in 1:9){ 
    for(enviro in c(1,3,6,9,12)){
      for(pool in c(0.0001,0.001,0.01,0.1,1)){
        for (disper in c("near","radius","global")) {
          studypool[indrow,1]<-replic
          studypool[indrow,2]<-patchi
          studypool[indrow,3]<-enviro
          studypool[indrow,4]<-pool
          studypool[indrow,5]<-disper
          tranb<-read.csv(paste("rep",replic,"map",36,disper,"env",enviro,"pool",pool,".csv"))
          trandata<-tail(tranb,1)
          DN<-round(apply(trandata,2,mean))
          DN<-DN[-1]
          S<-sum(DN)
          studypool[indrow,6]<-S
          lostind<-which(DN!=0)
          studypool[indrow,7]<-length(lostind)
          D<-0
          for (i in lostind) {
            D<-D+(DN[i]/S)^2
          }
          D<-round(D,4)
          studypool[indrow,8]<-1-D
          indrow<-1+indrow
        }
      }
    }
  }
}
write.csv(studypool,paste("/home/zoushuai/fragnew.csv"),row.names=F)





a<-read.csv(file.choose())

which(rowSums(a)!=144)
