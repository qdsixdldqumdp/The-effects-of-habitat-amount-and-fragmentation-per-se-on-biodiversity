#读取ha数据，按照p小于0.05的规则挑选AIC最小的拟合，并作图
library(figpatch)
library(ggplot2)
library(ggpubr)
library(Rmisc)
library(RColorBrewer)
library(stringr)
library(grid)
library(magick)
library(basicTrendline)
options(scipen = 10)


fity<-function(x,y,modelname,fnum){
  mymodel<-trendline_summary(x,y,model=modelname,summary=FALSE)# (y=a*x^b+c).
  a<-round(mymodel$parameter$a,2)
  b<-round(mymodel$parameter$b,2)
  deter_c<-ifelse(length(mymodel$parameter)>2,T,F)
  if(deter_c){c<-round(mymodel$parameter$c,2)}
  
  switch(fnum,
         fitting <- lm(y ~ x),#formula = "y = a*x + b"
         fitting <- lm(y ~ I(x^2) + x),#formula = "y = a*x^2 + b*x + c"
         fitting <- nls(y ~ SSexp2P(x, a, b), data = data.frame(x=x,y=y)),#formula = "y = a*exp(b*x)"
         fitting <- nls(y ~ SSexp3P(x, a, b, c), data = data.frame(x=x,y=y)), #formula = "y = a*exp(b*x) + c"
         fitting <- nls(y ~ SSpower3P(x, a, b, c),data = data.frame(x=x,y=y))#formula = "y = a*x^b + c"
         
  )
  
  predict(fitting)
}
pplotstandize<-function(p1){
  p1<-as.numeric(p1)
  if(p1<0.0001)
  {
    p2<-"~italic(p)~'<'~0.0001"
  }else if (p1<0.001)
  {
    p2<-"~italic(p)~'<'~0.001"
  }else if (p1<0.01)
  {
    p2<-"~italic(p)~'<'~0.01"
  }else if (p1<0.05)
  {
    p2<-"~italic(p)~'<'~0.05"
  }else{
    p1<-round(p1,4)
    ita<-as.character("~italic(p)~'='~")
    lin<-paste0(ita,p1)
    p2<-as.character(as.expression(lin))
  }
}
lunwen<-function(){
  mylunwen1<-"在测试divesity多样性指数对取样区域所在斑块面积变化的响应时设置了三个不同的扩散距离(D)，五个不同的迁入率(a)和五个不同的空间自相关程度，因此总共有75种组合。结果(下图)显示取样区域所在斑块大小与物种多样性无显著相关性的占比pnoef%（nnoef种参数组合，p ≥ 0.05）正相关占比为pposi% (nposi种参数组合，p<0.05)，负相关占比为pneg% (nneg种参数组合，p<0.05)，单峰曲线占比为pgao%（ngao种参数组合）。这表明大多数情况是支持栖息地总量假说的。"
  choosedata<- na.omit(modeldf)
  nsum<-length(choosedata$pvalue)
  if(diver=="sp"){
    diversity<-"物种丰富度"
  }else if(diver=="shannon"){
    diversity<-"香农威纳"
  }else if(diver=="simpson"){
    diversity<-"辛普森"
  }else if(diver=="invsimpson"){
    diversity<-"逆辛普森"
  }
  write.csv(choosedata,paste0("纯破碎化", diversity,".csv"))
  mylunwen1<-gsub("divesity",diversity,mylunwen1)
  noefdf<-subset(choosedata,pvalue>=0.05)
  length(noefdf$pvalue)->nnoef
  mylunwen1<-gsub("pnoef",round(nnoef/nsum*100,2),mylunwen1)
  mylunwen1<-gsub("nnoef",nnoef,mylunwen1)
  efdf<-subset(choosedata,pvalue<0.05)
  efdfline<-efdf[efdf$plotchoose=="choose1",]
  text<-efdfline$formula
  nneg<-sum( substr(text, 5, 5)== "-") 
  nposi<-length(text)-nneg
  mylunwen1<-gsub("pposi",round(nposi/nsum*100,2),mylunwen1)
  mylunwen1<-gsub("nposi",nposi,mylunwen1)
  mylunwen1<-gsub("pneg",round(nneg/nsum*100,2),mylunwen1)
  mylunwen1<-gsub("nneg",nneg,mylunwen1)
  efdfgao<-efdf[efdf$plotchoose!="choose1",]
  ngao<-length(efdfgao$pvalue)
  mylunwen1<-gsub("pgao",round(ngao/nsum*100,2),mylunwen1)
  mylunwen1<-gsub("ngao",ngao,mylunwen1)
  return(mylunwen1)
}





setwd("/Users/jingwufeng/Downloads/frag_test")
#file_name<-file.choose()
file_name<-"/Users/jingwufeng/Downloads/frag_test/patchmean0918.csv"
profragdata<-read.csv(file_name)#patchsize
profragdata$patch<-profragdata$patch+144
datarowname<-colnames(profragdata)
lendatarow<-length(datarowname)
stdatarowname<-which("reptime"==datarowname)
profragdata<-profragdata[,stdatarowname:lendatarow]

lunwenlist<-"a";lunwenind<-1





for (diver in c("sp","shannon","simpson" ,"invsimpson" )) {
  if(diver=="sp"){
  spautotitle<-"Species richness"
  }else if(diver=="shannon"){
    spautotitle<-"Shannon diversity"
  }else if(diver=="simpson"){
    spautotitle<-"Simpson diversity"
  }else if(diver=="invsimpson"){
    spautotitle<-"Invsimpson diversity"
  }
  
  letind<-1;myplotlist<-list()
  arrow <- arrow(length = unit(0.4, "cm"), type = "closed") 
  modelrow<-0
  modelvec<-c("line2P","line3P","exp2P","exp3P","power3P")
  modeldf<-data.frame(dispersal=1,spatial_auto=1,immigration=1,formula=1,pvalue=1,AIC=1,plotchoose=1)
  
  for (disper in c("near","radius","global")) {
    for(enviro in c(1,3,6,9,12)){
      envind<-which(enviro==profragdata$enviro)
      disind<-which(disper==profragdata$dispeisal)
      resultind<-intersect(disind,envind)
      data1<-profragdata[resultind,c(1,2,4,7)]
      attach(profragdata)
      data1$sp<-get(diver)[resultind]

      detach(profragdata)
      new_data <- summarySE(data1,measurevar = "sp",groupvars = c('patch','poolposible'))
      new_data$poolposible <- as.factor(new_data$poolposible)
      groupcol<-c(0.0001,0.001,0.01,0.1,1)
      p2<-NULL
      yposi<-NULL
      
      new_data2 <-new_data[,1:4] 
      
      ylim_max<-max(new_data$sp)
      ylim_min<-min(new_data$sp)
      for (colind in 1:5) {
        pool<-groupcol[colind]
        poolind<-which(pool==new_data$poolposible)
        yposi<-c(yposi,mean(new_data[poolind,4]))
        x<-new_data[poolind,1]
        y<-new_data[poolind,4]
        
        
        mymodel<-try(trendline_summary(x,y,model="line2P" ,summary = F))#(formula as: y=a*x+b)
        errordeter<-'try-error' %in% class(mymodel)
        if(errordeter) {
          modelrow<-1+modelrow
          modeldf[modelrow,1]<-disper
          modeldf[modelrow,2]<-enviro
          modeldf[modelrow,3]<-pool
          modeldf[modelrow,4]<-"y=a*x+b"
          modeldf[modelrow,5]<-1
          modeldf[modelrow,6]<-"aic"}else{
            a<-round(mymodel$parameter$a,2)
            b<-round(mymodel$parameter$b,2)
            myformula<-str_glue("y = {a}*x + {b}")
            pvalue<-mymodel$p.value
            aic<-mymodel$AIC
            modelrow<-1+modelrow
            modeldf[modelrow,1]<-disper
            modeldf[modelrow,2]<-enviro
            modeldf[modelrow,3]<-pool
            modeldf[modelrow,4]<-myformula
            modeldf[modelrow,5]<-pvalue
            modeldf[modelrow,6]<-aic
          }
        rowst<-modelrow#筛选p与aic
        
        mymodel<-try(trendline_summary(x,y,model="line3P",summary = F ))#(y=a*x^2+b*x+c)
        errordeter<-'try-error' %in% class(mymodel)
        if(errordeter) {
          modelrow<-1+modelrow
          modeldf[modelrow,1]<-disper
          modeldf[modelrow,2]<-enviro
          modeldf[modelrow,3]<-pool
          modeldf[modelrow,4]<-"y=a*x^2+b*x+c"
          modeldf[modelrow,5]<-1
          modeldf[modelrow,6]<-"aic"}else{
            a<-round(mymodel$parameter$a,2)
            b<-round(mymodel$parameter$b,2)
            c<-round(mymodel$parameter$c,2)
            myformula<-str_glue("y={a}*x^2+{b}*x+{c}")
            pvalue<-if(a==0){1}else{mymodel$p.value}
            aic<-mymodel$AIC
            modelrow<-1+modelrow
            modeldf[modelrow,1]<-disper
            modeldf[modelrow,2]<-enviro
            modeldf[modelrow,3]<-pool
            modeldf[modelrow,4]<-myformula
            modeldf[modelrow,5]<-pvalue
            modeldf[modelrow,6]<-aic}
        
        mymodel<-try(trendline_summary(x,y,model="exp2P",summary = F))# (y=a*exp(b*x))
        errordeter<-'try-error' %in% class(mymodel)
        if(errordeter) {
          modelrow<-1+modelrow
          modeldf[modelrow,1]<-disper
          modeldf[modelrow,2]<-enviro
          modeldf[modelrow,3]<-pool
          modeldf[modelrow,4]<-"y=a*exp(b*x)"
          modeldf[modelrow,5]<-1
          modeldf[modelrow,6]<-"aic"}else{
            a<-round(mymodel$parameter$a,2)
            b<-round(mymodel$parameter$b,2)
            myformula<-str_glue("y={a}*exp({b}*x)")
            pvalue<-if(a==0){1}else{mymodel$p.value}
            aic<-mymodel$AIC
            modelrow<-1+modelrow
            modeldf[modelrow,1]<-disper
            modeldf[modelrow,2]<-enviro
            modeldf[modelrow,3]<-pool
            modeldf[modelrow,4]<-myformula
            modeldf[modelrow,5]<-pvalue
            modeldf[modelrow,6]<-aic}
        
        mymodel<-try(trendline_summary(x,y,model="exp3P",summary = F))# (y=a*exp(b*x)+c)
        errordeter<-'try-error' %in% class(mymodel)
        if(errordeter) {
          modelrow<-1+modelrow
          modeldf[modelrow,1]<-disper
          modeldf[modelrow,2]<-enviro
          modeldf[modelrow,3]<-pool
          modeldf[modelrow,4]<-"a*exp(b*x)+c"
          modeldf[modelrow,5]<-1
          modeldf[modelrow,6]<-"aic"}else{
            a<-round(mymodel$parameter$a,2)
            b<-round(mymodel$parameter$b,2)
            c<-round(mymodel$parameter$c,2)
            myformula<-str_glue("y={a}*exp({b}*x)+{c}")
            pvalue<-if(a==0){1}else{mymodel$p.value}
            aic<-mymodel$AIC
            modelrow<-1+modelrow
            modeldf[modelrow,1]<-disper
            modeldf[modelrow,2]<-enviro
            modeldf[modelrow,3]<-pool
            modeldf[modelrow,4]<-myformula
            modeldf[modelrow,5]<-pvalue
            modeldf[modelrow,6]<-aic}
        
        
        rowen<-modelrow#筛选p与aic
        deterdf<-modeldf[rowst:rowen,4:6]
        deterdf[,1]<-1:4
        
          fnum<-1
          p2<-c(p2,deterdf[1,2])
          plotmodel<-modelvec[fnum]
          new_data2[poolind,4]<-fity(x,y,plotmodel,fnum)
          plot_row<-rowst-1+fnum
          modeldf[plot_row,7]<-"choose1"
          }
        
        
      
      length(p2)->p2_n
      for(p in 1:p2_n){
        p2[p]<-pplotstandize(p2[p])
      }
      
      
      col1 <- brewer.pal(5,"Dark2")
      attach(new_data)
      
      yfanwei<-ylim_max-ylim_min
      kexue<-format(signif(yfanwei,1),scientific = TRUE)
      zhengfu<-substr(kexue,3,3)
      beishu<- as.numeric(substr(kexue,4,5))
      beilv<-ifelse(zhengfu=="+",10^beishu,10^(-beishu))
      ifelse(yfanwei/beilv<5,kedu<-beilv,kedu<-beilv*2)
      myymin<-floor(ylim_min/beilv)*beilv
      
      yposilen<-ylim_max-myymin
      if(diver=="shannon"){if(disper=="near"){
        if(enviro==1){yposi<-c(3,13,seq(from=65,by=8,length.out=3))/100*yposilen+min(myymin)}else{yposi<-c(3,13,seq(from=55,by=8,length.out=3))/100*yposilen+min(myymin)}
        }else{
        yposi<-c(3,seq(from=55,by=8,length.out=4))/100*yposilen+min(myymin)}
      }else if(diver=="simpson"){if(disper=="near"){
        if(enviro<6){yposi<-c(3,11,19,75,85)/100*yposilen+min(myymin)
        }else if(enviro==12){yposi<-c(3,11,19,27,85)/100*yposilen+min(myymin)}else{yposi<-c(3,13,seq(from=75,by=8,length.out=3))/100*yposilen+min(myymin)}
      }else{ if(disper=="global"&enviro>1){yposi<-c(3,seq(from=62,by=8,length.out=4))/100*yposilen+min(myymin)}else{
        yposi<-c(3,13,seq(from=70,by=8,length.out=3))/100*yposilen+min(myymin)}
        }
      }else if(diver=="sp"){
        yposi<-seq(from=40,by=8,length.out=5)/100*yposilen+min(myymin)
      }else{
        yposi<-seq(from=30,by=8,length.out=5)/100*yposilen+min(myymin)
      }
      if(letind<1){
        posi<-c(0.25,0.6)}else{posi<-"none"}
      xzhiyu<-max(patch)-min(patch)
      plotsize<-10
      p_vx<-if(diver!="simpson"){rep(xzhiyu*0.7,each=5)}else{
        if(disper=="near"){c(xzhiyu*0.5,xzhiyu*0.5,rep(xzhiyu*0.7,each=3))}else{rep(xzhiyu*0.7,each=5)}}
      p_vx<-p_vx+min(patch)
      mybreak<-as.data.frame(table(profragdata$patch))[1,1]
      xcontinu<-0
      for (i in 1:length(mybreak)) {
        xcontinu[i]<-as.numeric(as.character(mybreak[i]))
      }
      
      mybreak<-as.data.frame(table(profragdata$patch))[,1]
      xcontinu<-as.numeric(as.character(mybreak))
      seqgap<-25
      xcontinu<-seq(round(min(xcontinu)/seqgap)*seqgap,max(xcontinu),seqgap)
      
 
      
      
      p<-
        ggplot(new_data,aes(x=patch,y=sp,group=poolposible,color=poolposible))+
        geom_point(size=plotsize/4) +
        scale_color_brewer(palette="Dark2")+
        geom_smooth(data = new_data,aes(x=patch,y=sp,group=poolposible,color=poolposible),size=plotsize/9,method = "lm", se=FALSE, formula = y~x)+
        geom_errorbar(size=plotsize/9,aes(ymin=sp-se, ymax=sp+se), width=xzhiyu*0.05)+
        annotate("text",x=p_vx,y=yposi, label = p2, parse = TRUE,size=plotsize,color=col1)+
        annotate("text",x=-Inf,y=Inf,hjust=-1,vjust=2,label=LETTERS[letind],color="black",fontface = "bold",size=plotsize)+
        xlab(NULL)+ylab(NULL)+#ylim(0,ylim_max)+
        guides(color = guide_legend(title = "Immigration rate"))+
        scale_x_continuous(breaks = xcontinu ,labels = if(letind%in%c(1:3*5)){paste0(xcontinu)}else{NULL})+
        scale_y_continuous( breaks =seq(myymin,ylim_max,kedu), labels =if(letind%in%c(1:5)){seq(myymin,ylim_max,kedu)}else{NULL},limits = c(myymin,ylim_max))+
        theme(
          panel.background = element_rect(fill = "white",color = "black",size=1.1),
          axis.text = element_text(size = 22,color="black",face="bold"),
          #axis.title.x = element_text(angle=0,hjust=1,vjust=1,size = 22,color="black",face="bold"),
          legend.key.height = unit(25, "pt"),
          legend.key.width = unit(25, "pt"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.position=posi) 
      
      
      detach(new_data)
      myplotlist[[letind]]<-p
      letind<-letind+1
      #print(LETTERS[letind])
    }}
  #write.csv(modeldf,"model-ps.csv",row.names = F)
  lunwenlist[lunwenind]<-lunwen()
  lunwenind<-lunwenind+1
  
  
  
  {
    jpgfile<-paste0("pho_v/size",diver,".png")
    png(jpgfile,width = 2400,height = 2400)
    
    sheji<-F
    title<-"The size of patch in which the sample area located"
    xlabtitle<-"Surrounding sample area piexls"
    layweight<-0.75
    
    
    #---------------------------------------------------
    grid.newpage()#paiban
    layweight<-0.75
    #x = pinx, y = piny, w = pinw, h = pinh, just = c("right", "bottom")
    pinx=1
    lengend_h<-0.02
    piny=0.04
    pinw=layweight*0.6
    pinh=layweight
    titlefontsize<-40
    textfontsize<-38
    auto_h<-0.03
    
    lengendpic<-image_read("/Users/jingwufeng/Downloads/frag_test/lengend.png")#图例--------------"red"
    lengendport<-viewport(x = pinx, y = piny, w = pinw, h = lengend_h,  just = c("right", "top"))
    grid.raster( lengendpic,vp=lengendport)
    #pushViewport(name = lengendport) 
    #if(sheji){grid.rect(gp = gpar(col = "green"))}else{print(lengendpic,vp = viewport(layout.pos.row = 1, layout.pos.col = 1))}
    #upViewport(1)
    
    #title
    ytitle<-piny
    vp1 <- viewport(x = pinx, y = ytitle-lengend_h, w = pinw, h = piny-lengend_h, just = c("right", "top"))  
    pushViewport(vp1) 
    if(sheji){grid.rect(gp = gpar(col = "red"))}
    if(!sheji){grid.text(title, x = 0.5, y = 0.7, gp = gpar(col = "black", fontfamily = "serif", fontsize = titlefontsize,fontface = "bold")) }
    upViewport(1)
    
    
    
    #Autocorrelation
    autoxposi<-pinx-pinw
    spx<-piny;spj<-layweight/5
    autovec<-c(1,3,6,9,12);autovec<-autovec[5:1]
    for (auto in 0:4) {
      vpau<-viewport(x = autoxposi, y = spx+spj*auto, w = auto_h, h =spj , just = c("right", "bottom")) 
      pushViewport(name = vpau)  
      if(!sheji){grid.text(paste0("A = ",autovec[auto+1]), x = 0.7, y = 0.5, rot=90, gp = gpar(col = "black", fontfamily = "serif", fontsize = textfontsize,fontface = "bold"))}
      if(sheji){grid.rect(gp = gpar(col = "blue"))}
      upViewport(1)
    }
    autotitlex<-autoxposi-auto_h
    vp2 <- viewport(x = autotitlex, y = spx, w = auto_h, h =layweight , just = c("right", "bottom")) 
    pushViewport(vp2)  
    #"Spatial autocorrelation"
    if(sheji){grid.rect(gp = gpar(col = "black"))}
    if(!sheji){grid.text(spautotitle, rot=90, x = 1, y = 0.5, gp = gpar(col = "black", fontfamily = "serif", fontsize = titlefontsize,fontface = "bold"))}
    upViewport(1)
    
    
    #dispersal
    spy<-piny+layweight;spx<-layweight/5
    autovec<-c("Global","Radius-5","Nearest-4")
    for (auto in 0:2) {
      v1<-viewport(x = 1-spx*auto, y = spy, w = spx, h = auto_h, just = c("right", "bottom"))
      pushViewport(name = v1)  
      if(sheji){grid.rect(gp = gpar(col = "green"))}
      if(!sheji){grid.text(paste0(autovec[auto+1]), x = 0.5, y = 0.3,gp = gpar(col = "black", fontfamily = "serif", fontsize = textfontsize,fontface = "bold")) }
      #grid.rect(gp = gpar(col = "green"))
      upViewport(1) 
    }
    v4<-viewport(x = pinx, y = spy+auto_h, w = layweight*0.6, h =auto_h , just = c("right", "bottom"))
    pushViewport(name = v4)  
    if(sheji){grid.rect(gp = gpar(col = "black"))}
    if(!sheji){grid.text("Dispersal distance", x = 0.5, y = 0,gp = gpar(col = "black", fontfamily = "serif", fontsize = titlefontsize,fontface = "bold")) }
    upViewport(1)
    
    
    
    
    vp<-viewport(x = pinx, y = piny, w = pinw, h = pinh, just = c("right", "bottom"))
    pushViewport(vp)
    if(sheji){grid.rect(gp = gpar(col = "yellow"))}
    if(!sheji){
      layout_1 <- grid.layout(nrow = 5, ncol = 3) 
      pushViewport(viewport(layout = layout_1)) 
      
      
      for (i in 1:3) {
        for (j in 1:5) {
          n<-i*5+j-5
          print(myplotlist[[n]], vp = viewport(layout.pos.row = j, layout.pos.col = i))
        }
      }
    }
    
    dev.off()
    #image_crop(image_read(jpgfile),geometry ="1300x2200+1200+420")
    image_write(image_crop(image_read(jpgfile),geometry ="1300x2200+1220+410"),jpgfile)#看图像大小
    
  }
  
}

jieguo<-data.frame("wu"=1,"hao"=1,"fu"=1,"qu"=1)
for (index in 1:4) {
  tiqu<-lunwenlist[index]
  weizhi<-str_locate_all(tiqu,"参数组合")|>as.data.frame()
  weizhi<-weizhi$start
  for (i in 1:length(weizhi)) {
    jieguo[index,i]<-str_extract_all(substr(tiqu,weizhi[i]-4,weizhi[i]-1), "\\d+")[[1]]|>as.numeric()
  }
}
zongjie<-"总体而言，在总共300中参数（75*4种多样性指数）组合中，取样区域所在斑块大小与多样性之间无相关的占比pnoef%（nnoef种参数组合，p ≥ 0.05）正相关占比为pposi% (nposi种参数组合，p<0.05)，负相关占比为pneg% (nneg种参数组合，p<0.05)，单峰曲线占比为pgao%（ngao种参数组合）。"
res_sum<-colSums(jieguo)|>as.vector()
tivec<-c("pnoef","nnoef","pposi","nposi","pneg","nneg","pgao","ngao")
for (ti in 1:4) {
  zongjie<-gsub(tivec[ti*2-1],res_sum[ti]/sum(res_sum)*100,zongjie)
  zongjie<-gsub(tivec[ti*2],res_sum[ti],zongjie)
}
lunwenlist
zongjie

