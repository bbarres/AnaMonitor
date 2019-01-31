###############################################################################
###############################################################################
#Different code for producing monitoring-related maps
###############################################################################
###############################################################################

#loading the packages necessary for the analysis
library(rgdal)
library(rgeos)
library(plotrix)
library(mapplots)
library(RColorBrewer)
library(dplyr)


###############################################################################
#loading the data
###############################################################################

#load geographical data of departements and regions
load("data/departeLight.RData")
load("data/regionsLight.RData")

#load the barycentre coordinates of departements and regions
load("data/coorddep.RData")
load("data/coordreg.RData")

#loading the data exported from PROSPER "tableau de bord"
datasampl<-read.table("data/2018_toutheme_cor.txt",header=TRUE,sep="\t")
datasampl<-merge(datasampl,coordreg,by.x="Region",by.y="reg_CODENAME")

#loading the result data by departement
dataresult<-read.table("data/2017_themefin.txt",header=TRUE,sep="\t")
#to streamline subsequent analysis, we turned the resistance status factor
#into two different columns
dataresult$rslt_RS<-factor(dataresult$rslt_RS,levels=c("R","S"))
dataresult$Resistant<-as.numeric(dataresult$rslt_RS==levels(dataresult$rslt_RS)[1])
dataresult$Sensitive<-as.numeric(dataresult$rslt_RS==levels(dataresult$rslt_RS)[2])
dataresult$Total<-dataresult$Resistant+dataresult$Sensitive


###############################################################################
#Monitoring sampling maps by regions
###############################################################################

#a note for futur-me: this code will probably need some minor updates in order
#to work properly on the raw exported-from-PROSPER file

#this script produce png images of sampling maps for each "Programme" listed 
#in the "datasampl" file
colovec<-brewer.pal(8,"Dark2")
for (i in 1:length(levels(datasampl$Programme))){
  temp<-datasampl[datasampl$Programme==levels(datasampl$Programme)[i],]
  png(file=paste("output/",temp$Bioagresseur,temp$Hote,temp$Pesticides,".png",
                 sep=""),width=8,height=4,units="in",res=300)
  op<-par(mfrow=c(1,2),mar=c(0,0,0,0))
  plot(regionsLight,lwd=3)
  title(main="Prélèvements attendus",line=-1)
  draw.pie(x=temp$longitude,y=temp$latitude,
           z=cbind((temp$Prel_attend),0),
           col=colovec[1],lty=0,
           radius=(sqrt(temp$Prel_attend)*30000),
           labels=NA)
  text(x=temp$longitude,y=temp$latitude,
       labels=as.character(temp$Prel_attend),cex=2)
  plot(regionsLight,lwd=3)
  title(main="Prélèvement(s) reçu(s)",line=-1)
  if (sum(temp$Prel_recep)==0) {
    points(x=temp$longitude,y=temp$latitude)
    text(x=temp$longitude,y=temp$latitude,
         labels=as.character(temp$Prel_recep),cex=2)
  } else {
    draw.pie(x=temp$longitude,y=temp$latitude,
             z=cbind((temp$Prel_recep),0),
             col=colovec[2],lty=0,
             radius=(sqrt(temp$Prel_recep)*30000),
             labels=NA)
    text(x=temp$longitude,y=temp$latitude,
         labels=as.character(temp$Prel_recep),cex=2)
  }
  par(op)
  dev.off()
}


###############################################################################
#Maps of the results of the monitoring
###############################################################################

#grouping and counting the resistant and sensitive samples
dataCamem<-dataresult %>% 
  group_by(themat_ID,dptmt,Resistant) %>% 
  summarise(Resist=sum(Resistant),Sensi=sum(Sensitive),Tot=sum(Total))

data2map<-merge(dataCamen,coorddep,by.x="dptmt",by.y="dep_ID")

op<-par(mar=c(0,0,0,0))
plot(departeLight,col="grey70")
draw.pie(x=data2map$longitude,y=data2map$latitude,
         z=cbind((as.numeric(as.character(data2map$Resist))),
                 (as.numeric(as.character(data2map$Sensi)))),
         col=c("red","blue"),
         radius=(sqrt(as.numeric(as.character(data2map$Tot)))*30000),
         labels=NA)
par(op)



###############################################################################
#END
###############################################################################