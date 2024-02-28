library(stringr)
setwd("/home/zoushuai/patchsize0708")
formtxt<-readLines("form.R")
for(i in c(112,12,128,25,38,52,66,81,96,0)){
  dir.create(str_glue("/home/zoushuai/patchsize0708/{i}"))
  formtxt[1]<-paste("setwd(\"",str_glue("/home/zoushuai/patchsize0708/{i}"),"\")",sep = "")
  formtxt[160]<-str_glue("path<-{i}")
  formtxt[30]<-paste0("source(\"","/home/zoushuai/patchsize0708/disthemap/",str_glue("map{i}.R"),"\")")
  formtxt[162]<-"for(replic in 1:25){"
  write.table(formtxt,str_glue("size{i}.R"),row.names = F,quote = F,col.names = F)
  formtxt[162]<-"for(replic in 25:50){"
  write.table(formtxt,str_glue("size{i}_1.R"),row.names = F,quote = F,col.names = F)
 
  
}
  
    
    
    
sh<-"a";shind<-1   
for(i in c(112,12,128,25,38,52,66,81,96,0)){    
sh[shind]<-str_glue("nohup R  -f /home/zoushuai/patchsize0708/size{i}.R&")
sh[shind+10]<-str_glue("nohup R  -f /home/zoushuai/patchsize0708/size{i}_1.R&")
shind<-shind+1


}
write.table(sh,"patchsize.sh",row.names = F,quote = F,col.names = F)
    