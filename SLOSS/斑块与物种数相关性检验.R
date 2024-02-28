setwd("/Users/jingwufeng/Downloads/habitat fragment/SLOSS/fragtest")
library(dplyr)
#abundata<-read.csv(file.choose())
#metadata<-read.csv(file.choose())
abundata<-read.csv("/Users/jingwufeng/Downloads/habitat fragment/SLOSS/doi_10.5061_dryad.595718c__v1/abundances_utf8.csv")
metadata<-read.csv( "/Users/jingwufeng/Downloads/habitat fragment/SLOSS/doi_10.5061_dryad.595718c__v1/metadata_utf8.csv")
meta<-metadata[,c("refshort","site_id" ,"site_size" )]
meta<-meta[!meta$site_size=="continuous",]
meta$site_size<-as.numeric(meta$site_size)
meta<-aggregate(site_size~refshort+ site_id,meta,mean)
wenxian<-as.character(as.data.frame(table(abundata$refshort))$Var1)
xiangguan<-data.frame()
for (i in 1:length(wenxian)) {
  wenxian_frame<-abundata[abundata$refshort==wenxian[i],]
  site_sp<-data.frame()
  plot_list<-as.character(as.data.frame(table(wenxian_frame$site_id))$Var1)
  for (j in 1:length(plot_list)) {
    n_sp<-wenxian_frame[wenxian_frame$site_id==plot_list[j],]
    site_sp[j,1]<-plot_list[j]
    site_sp[j,2]<-length(table(n_sp$scientific_name))
    ind1<-meta[which(meta$refshort==wenxian[i]),] 
    if(all(!ind1$site_id==plot_list[j])){site_sp[j,3]<-max(ind1$site_size)*2}else{
    site_sp[j,3]<-ind1[ind1$site_id==plot_list[j],3]}
  }
  colnames(site_sp)<-c("plot_id","sp","size")
  write.csv(site_sp,paste0(wenxian[i],".csv"))
  xiangguan[i,1]<-wenxian[i]
  if(length(site_sp$sp)>3){
  #xiangguan[i,2]<-cor(site_sp[,-1], method="spearman")[2]
  xiangguan[i,2]<-cor(site_sp[,-1], method="pearson")[2]
  }else{xiangguan[i,2]<-NA}
}

write.csv(xiangguan,"/Users/jingwufeng/Downloads/habitat fragment/SLOSS/pearson.csv")

zongjie<-na.omit(xiangguan)
zongjie[abs(zongjie$V2)>=0.9,]->xianzhu
zongjie[abs(zongjie$V2)<0.9,]->buxianzhu
length(xianzhu$V1);length(buxianzhu$V1)
length(xianzhu$V1)/length(zongjie$V1)
length(buxianzhu$V1)/length(zongjie$V1)

zongjie[abs(zongjie$V2)>=0.7,]->xianzhu
zongjie[abs(zongjie$V2)<0.7,]->buxianzhu
length(xianzhu$V1);length(buxianzhu$V1);length(zongjie$V1)
length(xianzhu$V1)/length(zongjie$V1)
length(buxianzhu$V1)/length(zongjie$V1)

