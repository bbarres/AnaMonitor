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

#load geographical data of regions
load("data/REG_SHP.RData")
#load the barycentre coordinates of regions
load("data/coordreg.RData")

#loading the data exported from PROSPER "tableau de bord"
sampMoni<-read.table("data/testdata.txt",header=TRUE,
                     sep="\t",stringsAsFactors=TRUE)
sampMoni<-merge(sampMoni,coordreg,by.x="Region",by.y="reg_CODENAME")

#defining a set of colors
colovec<-brewer.pal(8,"Dark2")


##############################################################################/
#defining additional function for the mapping####
##############################################################################/

#function for a scale, found in "Auxiliary Cartographic Functions in R: 
#North Arrow, Scale Bar, and Label with a Leader Arrow", Tanimura et al 2007, 
#J of Statistical software
#The code has been slightly modified in order to convert the meter in km
scalebar <- function(loc,length,unit="km",division.cex=.8,...) {
  if(missing(loc)) stop("loc is missing")
  if(missing(length)) stop("length is missing")
  x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1]
  y <- c(0,length/(10*3:1))+loc[2]
  cols <- rep(c("black","white"),2)
  for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
  for (i in 1:5) segments(x[i],y[2],x[i],y[3])
  labels <- (x[c(1,3)]-loc[1])/1000
  labels <- append(labels,paste((x[5]-loc[1])/1000,unit))
  text(x[c(1,3,5)],y[4],labels=labels,adj=c(0.5,0),cex=division.cex)
}


##############################################################################/
#Mapping planned and realized sampling####
##############################################################################/

#this script produce png images of sampling maps
temp<-sampMoni
png(file=paste("output/",temp$Bioagresseur,temp$Hote,temp$Pesticides,".png",
               sep=""),width=8,height=4,units="in",res=300)
op<-par(mfrow=c(1,2),mar=c(0,0,0,0))
plot(REG_SHP,lwd=3)
title(main="Prélèvement(s) attendu(s)",line=-1)
draw.pie(x=temp$longitude,y=temp$latitude,
         z=cbind((temp$Prel_attend),0),
         col=colovec[1],lty=0,
         radius=(sqrt(temp$Prel_attend)*30000),
         labels=NA)
text(x=temp$longitude,y=temp$latitude,
     labels=as.character(temp$Prel_attend),cex=2)
scalebar(c(191260,6060000),300000,"km",division.cex=0.8)
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