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
par(family = "serif")


guass_peak<-function(x,y){
  eDigit<-5
  my_model <- function(x, y0, a, x0, b) {  
    y0 + a * exp(-0.5 * ((x - x0) / b)^2)  
  }  
  
  n = length(x)
  k = 2
  
  initial_y0 <- mean(y) # 初始y0值设为y的均值  
  initial_a <- max(y) - initial_y0 # 初始a值设为y的最大值与y0的差  
  initial_x0 <- mean(x) # 初始x0值设为x的均值  
  initial_b <- sd(x) # 初始b值设为x的标准差  
  #mymodel<-try(trendline_summary(x,y,model="line3P",summary = F ))#(y=a*x^2+b*x+c)
  #errordeter<-'try-error' %in% class(mymodel)
  
  fit <- try(nls2::nls2(y ~ my_model(x, y0, a, x0, b),   
                        start = list(y0 = initial_y0, a = initial_a, x0 = initial_x0, b = initial_b),   
                        algorithm = "port"  ) )
  errordeter<-'try-error' %in% class(fit)
  if(!errordeter){
    sum.exp2P <- summary(fit)
    ss.res <- sum((residuals(fit))^2)
    ss.total.uncor <- sum(y^2)
    ss.total.cor <- sum((y - mean(y))^2)
    ss.reg <- ss.total.cor - ss.res
    dfR = k - 1
    dfE = n - k
    Fval = (ss.reg/dfR)/(ss.res/dfE)
    pval = pf(Fval, dfR, dfE, lower.tail = F)
    pval <- unname(pval)
    #RSE <- sum.exp2P$sigma
    #SSE <- (RSE^2) * (n - 1)
    #adjr2 <- 1 - SSE/((var(y)) * (n - 1))
    #r2 <- 1 - (1 - adjr2) * ((n - k)/(n - 1))
    coeff = coef(sum.exp2P)
    paralist<-setNames(as.list(coeff[,1]), rownames(coeff))
    
    #cat("\nF-statistic:", format(Fval, digits = eDigit),  "on", dfR, "and", dfE, "DF, ", "p-value:", format(pval, digits = eDigit), "\n")
    AIC = as.numeric(format(AIC(fit), digits = eDigit))
    fit_result<-list(AIC=AIC,p=pval,parameter=paralist)
  }else{fit_result<-NA}
  
  return(fit_result)
}
fity<-function(x,y,modelname,fnum){
 
  if(modelname!="guass" ){
  mymodel<-trendline_summary(x,y,model=modelname,summary=FALSE)# (y=a*x^b+c).
  a<-round(mymodel$parameter$a,2)
  b<-round(mymodel$parameter$b,2)
  deter_c<-ifelse(length(mymodel$parameter)>2,T,F)
  if(deter_c){c<-round(mymodel$parameter$c,2)}
  
  switch(fnum,
         fitting <- lm(y ~ x),#formula = "y = a*x + b"
         #fitting <- lm(y ~ I(x^2) + x),#formula = "y = a*x^2 + b*x + c"
         fitting <- nls(y ~ SSexp2P(x, a, b), data = data.frame(x=x,y=y)),#formula = "y = a*exp(b*x)"
         fitting <- nls(y ~ SSexp3P(x, a, b, c), data = data.frame(x=x,y=y)), #formula = "y = a*exp(b*x) + c"
         #fitting <- nls(y ~ SSpower3P(x, a, b, c),data = data.frame(x=x,y=y)),#formula = "y = a*x^b + c"
         
  )}else{ 
    my_model <- function(x, y0, a, x0, b) {  
    y0 + a * exp(-0.5 * ((x - x0) / b)^2)  
  } 
    initial_y0 <- mean(y) # 初始y0值设为y的均值  
    initial_a <- max(y) - initial_y0 # 初始a值设为y的最大值与y0的差  
    initial_x0 <- mean(x) # 初始x0值设为x的均值  
    initial_b <- sd(x) # 初始b值设为x的标准差  
  fitting <-nls2::nls2(y ~ my_model(x, y0, a, x0, b),   
                       start = list(y0 = initial_y0, a = initial_a, x0 = initial_x0, b = initial_b),   
                       algorithm = "port"  )#formula = "y =y0 + a * exp(-0.5 * ((x - x0) / b)^2) "
  }
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
  mylunwen1<-"在测试divesity多样性指数对纯破碎化程度变化的响应时设置了三个不同的扩散距离(D)，五个不同的迁入率(a)和五个不同的空间自相关程度，因此总共有75种组合。结果显示取样区域所在斑块大小与物种多样性无显著相关性的占比pnoef%（nnoef种参数组合，p≥0.05），正相关占比为pposi% (nposi种参数组合，p<0.05)，负相关占比为pneg% (nneg种参数组合，p<0.05)，单峰曲线占比为pgao%（ngao种参数组合）。"
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
  #write.csv(choosedata,paste0("纯破碎化", diversity,".csv"))
  mylunwen1<-gsub("divesity",diversity,mylunwen1)
  noefdf<-subset(choosedata,pvalue>=0.05)
  length(noefdf$pvalue)->nnoef
  mylunwen1<-gsub("pnoef",round(nnoef/nsum*100,2),mylunwen1)
  mylunwen1<-gsub("nnoef",nnoef,mylunwen1)
  efdf<-subset(choosedata,pvalue<0.05)
  efdfline<-efdf[efdf$plotchoose=="choose1",]
  text<-efdfline$formula
  nneg<-sum( substr(text, 5, 5)== "-") 
  efdfmi<-efdf[efdf$plotchoose=="choose2"|efdf$plotchoose=="choose3",]
  nmi<-length(efdfmi$dispersal)
  nposi<-length(text)-nneg+nmi
  mylunwen1<-gsub("pposi",round(nposi/nsum*100,2),mylunwen1)
  mylunwen1<-gsub("nposi",nposi,mylunwen1)
  mylunwen1<-gsub("pneg",round(nneg/nsum*100,2),mylunwen1)
  mylunwen1<-gsub("nneg",nneg,mylunwen1)
  efdfgao<-efdf[efdf$plotchoose=="choose4",]
  ngao<-length(efdfgao$pvalue)
  mylunwen1<-gsub("pgao",round(ngao/nsum*100,2),mylunwen1)
  mylunwen1<-gsub("ngao",ngao,mylunwen1)
  return(mylunwen1)
}

#file_name<-file.choose()
file_name<-"/Users/jingwufeng/Downloads/frag_test/fragmean0918.csv"
profragdata<-read.csv(file_name)#frag
setwd("/Users/jingwufeng/Downloads/frag_test")
datarowname<-colnames(profragdata)
lendatarow<-length(datarowname)
stdatarowname<-which("reptime"==datarowname)
profragdata<-profragdata[,stdatarowname:lendatarow]
lunwenlist<-"a";lunwenind<-1



for (diver in c("sp","shannon","simpson" ,"invsimpson" )) {# diver="sp" 
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
  modelvec<-c("line2P","line3P","exp2P","exp3P","power3P","guass")[c(1,3,4,6)]
  modeldf<-data.frame(dispersal=1,spatial_auto=1,immigration=1,formula=1,pvalue=1,AIC=1,plotchoose=1)
  
  for (disper in c("near","radius","global")) {# disper="radius"
    for(enviro in c(1,3,6,9,12)){# enviro=9
      envind<-which(enviro==profragdata$enviro)
      disind<-which(disper==profragdata$dispeisal)
      resultind<-intersect(disind,envind)
      data1<-profragdata[resultind,c(1,2,4,7)]
      attach(profragdata)
      data1$sp<-get(diver)[resultind]
      ylim_max<-max(data1$sp)
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
            a<-round(mymodel$parameter$a,4)
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
            a<-round(mymodel$parameter$a,4)
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
            a<-round(mymodel$parameter$a,4)
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
        
        mymodel<-guass_peak(x,y)
        errordeter<-length(mymodel)==1
        if(errordeter) {
          modelrow<-1+modelrow
          modeldf[modelrow,1]<-disper
          modeldf[modelrow,2]<-enviro
          modeldf[modelrow,3]<-pool
          modeldf[modelrow,4]<-"y0 + a * exp(-0.5 * ((x - x0) / b)^2) "
          modeldf[modelrow,5]<-1
          modeldf[modelrow,6]<-"aic"
          }else{
            a<-round(mymodel$parameter$a,4)
            b<-round(mymodel$parameter$b,2)
            x0<-round(mymodel$parameter$x0,2)
            y0<-round(mymodel$parameter$y0,2)
            myformula<-str_glue("y={y0} + {a} * exp(-0.5 * ((x - {x0}) / {b})^2) ")
            pvalue<-if(a==0){1}else{mymodel$p}
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
        #if(disper=="near"){moxing<-deterdf[1,2]<0.05 }else{moxing<-T}
        moxing<-deterdf[1,2]<0.05
        if(moxing){
          fnum<-1
          p2<-c(p2,deterdf[1,2])
          plotmodel<-modelvec[fnum]
          new_data2[poolind,4]<-fity(x,y,plotmodel,fnum)
          plot_row<-rowst-1+fnum
          modeldf[plot_row,7]<-"choose1"}else{
            deterdf_p<-subset(deterdf,pvalue<0.05)
            if(length(deterdf_p$formula)==0){
              #print(rowen)
              fnum<-deterdf[which.min(deterdf$pvalue),1]
              p2<-c(p2,deterdf[which.min(deterdf$pvalue),2])
              plotmodel<-modelvec[fnum]
              new_data2[poolind,4]<-fity(x,y,plotmodel,fnum)
              plot_row<-rowst-1+fnum
              modeldf[plot_row,7]<-paste0("choose",fnum)
            }else{
              fnum<-deterdf_p[which.min(deterdf_p$AIC),1]
              p2<-c(p2,deterdf_p[which.min(deterdf_p$AIC),2])
              plotmodel<-modelvec[fnum]
              new_data2[poolind,4]<-fity(x,y,plotmodel,fnum)
              plot_row<-rowst-1+fnum
              modeldf[plot_row,7]<-paste0("choose",fnum)
            }
          }
        
        
        
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
      if(diver=="shannon"){
        if(disper=="near"){yposi<-seq(from=5,by=8,length.out=5)/100*yposilen+min(myymin)}else{
        yposi<-c(5,seq(from=65,by=8,length.out=4))/100*yposilen+min(myymin)}
      }else if(diver=="simpson"){
        if(disper=="near"){yposi<-seq(from=5,by=8,length.out=5)/100*yposilen+min(myymin)}else{
        yposi<-c(3,11,19,seq(from=85,by=8,length.out=2))/100*yposilen+min(myymin)}
      }else if(diver=="sp"){ 
        if(disper=="near"){yposi<-seq(from=50,by=8,length.out=5)/100*yposilen+min(myymin)}else{
        yposi<-seq(from=40,by=8,length.out=5)/100*yposilen+min(myymin)}
      }else{
        if(disper=="near"){yposi<-seq(from=50,by=8,length.out=5)/100*yposilen+min(myymin)}else{
        yposi<-seq(from=30,by=8,length.out=5)/100*yposilen+min(myymin)}
      }
      if(letind<1){
        posi<-c(0.25,0.6)}else{posi<-"none"}
      
      xzhiyu<-max(patch)-min(patch)
      p_vx<-if(diver=="simpson"){
        if(disper=="near"){rep(xzhiyu*0.7,each=5)}else if(disper=="radius"){rep(xzhiyu*0.4,each=5)}else{rep(xzhiyu*0.5,each=5)}
      }else if(diver=="shannon"){if(disper=="near"){rep(xzhiyu*0.8,each=5)}else{rep(xzhiyu*0.7,each=5)}
      }else if(diver=="sp"){rep(xzhiyu*0.35,each=5)
      }else{rep(xzhiyu*0.4,each=5)
      }
      
      plotsize<-10
      p<-
        ggplot(new_data,aes(x=patch,y=sp,group=poolposible,color=poolposible))+
        geom_point(size=plotsize/4) +
        scale_color_brewer(palette="Dark2")+
        geom_line(data = new_data2,aes(x=patch,y=sp,group=poolposible,color=poolposible),size=1)+
        geom_errorbar(size=plotsize/9,aes(ymin=sp-se, ymax=sp+se), width=xzhiyu*0.05)+
        annotate("text",x=p_vx,y=yposi, label = p2, parse = TRUE,size=plotsize,color=col1)+
        annotate("text",x=-Inf,y=Inf,hjust=-0.5,vjust=1.1,label=LETTERS[letind],color="black",fontface = "bold",size=plotsize)+
        xlab(NULL)+ylab(NULL)+#ylim(0,ylim_max)+
        guides(color = guide_legend(title = "Immigration rate"))+
       scale_x_continuous(breaks =c(0,20,40,60)  ,labels =if(letind%in%c(1:3*5)){ paste0(c(0,20,40,60))}else{NULL})+#c(1,2,4,8,16,36,72)
        scale_y_continuous( breaks =seq(myymin,ylim_max,kedu), labels =if(letind%in%c(1:5)){seq(myymin,ylim_max,kedu)}else{NULL},limits = c(myymin,ylim_max*1.04))+
        theme(
          panel.background = element_rect(fill = "white",color = "black",size=1.1),
          axis.text = element_text(size = 22,color="black",face="bold"),
          #axis.title.x = element_text(angle=0,hjust=1,vjust=1,size = 22,color="black",face="bold"),
          legend.key.height = unit(25, "pt"),
          legend.key.width = unit(25, "pt"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.position=posi) 
      #p
      
      detach(new_data)
      myplotlist[[letind]]<-p
      letind<-letind+1
      #print(LETTERS[letind])
    }}
  #write.csv(modeldf,paste0(diver,"model-frag.csv"),row.names = F)
  lunwenlist[lunwenind]<-lunwen()
  lunwenind<-lunwenind+1
  
  
  
  {  
    jpgfile<-paste0("pho_v/frag",diver,".png")
    png(jpgfile,width = 2400,height = 2400)
    
    sheji<-F
    title<-"Fragmentation per se (measured by the number of patches)"
    xlabtitle<-"Number of patches"
    
    
    
    
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
zongjie<-"总体而言，在总共300中参数（75*4种多样性指数）组合中，纯破碎化程度与多样性之间无相关的占比pnoef%（nnoef种参数组合，p ≥ 0.05）正相关占比为pposi% (nposi种参数组合，p<0.05)，负相关占比为pneg% (nneg种参数组合，p<0.05)，单峰曲线占比为pgao%（ngao种参数组合）。"
res_sum<-colSums(jieguo)|>as.vector()
tivec<-c("pnoef","nnoef","pposi","nposi","pneg","nneg","pgao","ngao")
for (ti in 1:4) {
  zongjie<-gsub(tivec[ti*2-1],res_sum[ti]/sum(res_sum)*100,zongjie)
  zongjie<-gsub(tivec[ti*2],res_sum[ti],zongjie)
}
lunwenlist
zongjie


