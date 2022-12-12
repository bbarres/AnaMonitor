##############################################################################/
##############################################################################/
#Script to produce the monitoring sampling
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(mapplots)
library(RColorBrewer)
library(sp)
library(dplyr)


##############################################################################/
#loading the data####
##############################################################################/

#load geographical data of departements and regions
load("data/DEP_SHP.RData")
load("data/REG_SHP.RData")

#load the barycentre coordinates of departements and regions
load("data/coorddep.RData")
load("data/coordreg.RData")

#loading the data exported from PROSPER "tableau de bord"
rezMoni<-read.table("data/2022_plasmo_rez.txt",header=TRUE,sep="\t",
                    colClasses="character")
#turning the resistance status factor into two different columns
rezMoni$rslt_RS<-factor(rezMoni$rslt_RS,levels=c("R","S"))
rezMoni$Resistant<-as.numeric(rezMoni$rslt_RS==levels(rezMoni$rslt_RS)[1])
rezMoni$Sensitive<-as.numeric(rezMoni$rslt_RS==levels(rezMoni$rslt_RS)[2])
rezMoni$Total<-rezMoni$Resistant+rezMoni$Sensitive
#grouping and counting the resistant and sensitive samples
rezMoni$themat_ID<-paste(rezMoni$themat_ID,rezMoni$SA)
rezMoni$themat_ID<-as.factor(rezMoni$themat_ID)
dataCamem<-rezMoni %>% 
  group_by(themat_ID,dptmt,pest,host) %>% 
  summarise(Resist=sum(Resistant),Sensi=sum(Sensitive),Tot=sum(Total))
data2map<-merge(dataCamem,coorddep,by.x="dptmt",by.y="dep_ID")

#defining a set of colors
colovec<-c(brewer.pal(9,"Reds")[7],brewer.pal(9,"Blues")[7])


##############################################################################/
#Maps of the results of the monitoring####
##############################################################################/

#mapping the results for each "programme"
colovec<-c(brewer.pal(9,"Reds")[7],brewer.pal(9,"Blues")[7])
for (i in 1:length(levels(data2map$themat_ID))){
  temp<-data2map[data2map$themat_ID==levels(data2map$themat_ID)[i],]
  png(file=paste("output/",temp$themat_ID,temp$pest,".png",sep=""),
      width=4,height=4,units="in",res=300)
  op<-par(mar=c(0,0,0,0))
  plot(DEP_SHP,border="grey70")
  plot(REG_SHP,lwd=2,add=TRUE)
  draw.pie(x=temp$longitude,y=temp$latitude,
           z=cbind((as.numeric(as.character(temp$Resist))),
                   (as.numeric(as.character(temp$Sensi)))),
           col=colovec,lty=0,
           radius=(sqrt(as.numeric(as.character(temp$Tot)))*15000),
           labels=NA)
  text(x=temp$longitude,y=temp$latitude,
       labels=as.character(temp$Tot),cex=1.2)
  #scalebar(c(191260,6060000),300000,"km",division.cex=0.8)
  par(op)
  dev.off()
}


##############################################################################/
#Distribution of the % of the ####
##############################################################################/

#plot histogram for the distribution of the percentage of sporulation
for (i in 1:length(levels(rezMoni$themat_ID))){
  temp<-rezMoni[rezMoni$themat_ID==levels(rezMoni$themat_ID)[i],]
  png(file=paste("output/","hist",temp$themat_ID,temp$pest,".png",sep=""),
      width=4,height=7,units="in",res=300)
  hist(as.numeric(temp$Perc_spor),ann=FALSE,las=1,freq=TRUE,
       include.lowest=TRUE,breaks=seq(0,100,by=10))
  dev.off()
}


##############################################################################/
#END
##############################################################################/