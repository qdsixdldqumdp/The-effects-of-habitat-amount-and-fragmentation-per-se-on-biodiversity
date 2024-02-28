library(stringr)
#read.table()函数会默认把双引号引起来的内容识别为一部分，使用quote = ""禁用该功能
#由于read.table()函数会默认把#后面的字符识别为注释，使用comment.char = ""禁用该功能。

setwd("/home/zoushuai/frag")

formtxt<-readLines("fragform.R")
patchivec<-c(1, 2, 4, 8, 16, 36, 72)

for(patchi in patchivec){
  #dir.create(str_glue("/home/zoushuai/frag/{patchi}"))
  formtxt[1]<-paste("setwd(\"",str_glue("/home/zoushuai/frag/{patchi}"),"\")",sep = "")
  formtxt[176]<-str_glue("patchi<-{patchi}")
  write.table(formtxt,str_glue("frag{patchi}.R"),row.names = F,quote = F,col.names = F)
}