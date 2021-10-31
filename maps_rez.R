##############################################################################/
##############################################################################/
#Script to produce the monitoring sampling
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(mapplots)
library(RColorBrewer)
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
rezMoni<-read.table("data/testdatarez.txt",header=TRUE,sep="\t",
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
#Maps of the results of the monitoring####
##############################################################################/

#mapping the results for each "programme" of the monitoring
  temp<-data2map
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
  scalebar(c(191260,6060000),300000,"km",division.cex=0.8)
  par(op)
  dev.off()


##############################################################################/
#END
##############################################################################/