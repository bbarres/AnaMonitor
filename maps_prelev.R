##############################################################################/
##############################################################################/
#Script to produce the monitoring sampling
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(rgdal)
library(rgeos)
library(plotrix)
library(mapplots)
library(RColorBrewer)
library(dplyr)


##############################################################################/
#loading the data####
##############################################################################/

#load geographical data of regions
load("data/REG_SHP.RData")
#load the barycentre coordinates of regions
load("data/coordreg.RData")

#loading the data exported from PROSPER "tableau de bord"
sampMoni<-read.table("data/testdata.txt",header=TRUE,
                     sep="\t",stringsAsFactors=TRUE)
sampMoni<-merge(sampMoni,coordreg,by.x="Region",
                by.y="reg_CODENAME")

#defining a set of colors
colovec<-brewer.pal(8,"Dark2")


##############################################################################/
#Mapping planned and realized sampling####
##############################################################################/

#this script produce png images of sampling maps for each "Programme" listed 
#in the "datasampl" file
temp<-sampMoni
png(file=paste("output/",temp$Bioagresseur,temp$Hote,temp$Pesticides,".png",
               sep=""),width=8,height=4,units="in",res=300)
op<-par(mfrow=c(1,2),mar=c(0,0,0,0))
plot(REG_SHP,lwd=3)
title(main="Prélèvements attendus",line=-1)
draw.pie(x=temp$longitude,y=temp$latitude,
         z=cbind((temp$Prel_attend),0),
         col=colovec[1],lty=0,
         radius=(sqrt(temp$Prel_attend)*30000),
         labels=NA)
text(x=temp$longitude,y=temp$latitude,
     labels=as.character(temp$Prel_attend),cex=2)
plot(REG_SHP,lwd=3)
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


##############################################################################/
#END
##############################################################################/