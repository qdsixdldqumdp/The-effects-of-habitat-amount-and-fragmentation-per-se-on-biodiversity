library(stringr)
#read.table()函数会默认把双引号引起来的内容识别为一部分，使用quote = ""禁用该功能
#由于read.table()函数会默认把#后面的字符识别为注释，使用comment.char = ""禁用该功能。

#dir.create("/home/zoushuai/ha")
setwd("/home/zoushuai/ha")

formtxt<-readLines("form.R")
hamapxy<-c(0, 52, 112, 180, 256, 340, 432, 532, 640, 756 , 880 )
hamapvec<-sqrt(hamapxy+144)

for(hamaplen in hamapvec){
  dir.create(str_glue("/home/zoushuai/ha/{hamaplen}"))
  formtxt[1]<-paste("setwd(\"",str_glue("/home/zoushuai/ha/{hamaplen}"),"\")",sep = "")
  formtxt[184]<-str_glue("inisimu(80,{hamaplen},{hamaplen},enviro)")
  write.table(formtxt,str_glue("ha{hamaplen}.R"),row.names = F,quote = F,col.names = F)
}














for(replicnumber in 1:50){
  formtxt[1]<-paste("setwd(\"",str_glue("/home/zoushuai/ha/{hamaplen}"),"\")",sep = "")
  formtxt[178]<-str_glue("inisimu(80,{hamaplen},{hamaplen},enviro)")
  formtxt[174]<-str_glue("replic<-{replicnumber}")
  write.table(formtxt,str_glue("ha{replicnumber}.R"),row.names = F,quote = F,col.names = F)
}
}









a<-read.table("form.txt",check.names=F,sep="\n",quote = "",comment.char = "");head(a)
a

write.table(lines,"qaq.R",row.names = F,quote = F,col.names = F)
lines <- readLines("form.txt")
lines[174]<-str_glue("replic<-{replicnumber}")

str_replace_all(lines[174],fixed("replic<-1"),str_glue("replic<-{replicnumber}"))




write.table(a,"qaq.R",row.names = F,quote = F,col.names = F)
replicnumber<-10
str_replace_all(a,fixed("replic<-1"),str_glue("replic<-{replicnumber}"))





str_replace_all("replic<-1 wj nasjnk",fixed("replic<-1"),str_glue("replic<-{replicnumber}"))
