##############################################################################/
##############################################################################/
#Code for the production of the map for the Rothamsted talk
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(rgdal)
library(rgeos)
library(plotrix)
library(mapplots)
library(RColorBrewer)
library(dplyr)
library(gdata)


##############################################################################/
#loading the data####
##############################################################################/

#load geographical data of departements and regions
load("data/departeLight.RData")
load("data/regionsLight.RData")

#load the barycentre coordinates of departements and regions
load("data/coorddep.RData")
load("data/coordreg.RData")

#loading the result data by departement
dataresult<-read.table("data/mildiou_Rotham.txt",header=TRUE,sep="\t")
#to streamline subsequent analysis, we turned the resistance status factor
#into two different columns
dataresult$rslt_RS<-factor(dataresult$rslt_RS,levels=c("R","S","RS"))
dataresult$Resistant<-as.numeric(dataresult$rslt_RS==
                                   levels(dataresult$rslt_RS)[1])
dataresult$Sensitive<-as.numeric(dataresult$rslt_RS==
                                   levels(dataresult$rslt_RS)[2])
dataresult$Suspect<-as.numeric(dataresult$rslt_RS==
                                 levels(dataresult$rslt_RS)[3])
dataresult$Total<-dataresult$Resistant+dataresult$Sensitive+dataresult$Suspect



##############################################################################/
#Maps of the results of the monitoring####
##############################################################################/

#a note for futur-me: this code will probably need some minor updates in order
#to work properly on the raw exported-from-PROSPER file

#grouping and counting the resistant and sensitive samples
dataCamem<-dataresult[dataresult$year==2012,] %>% 
  group_by(themat_ID,dptmt,pest,host) %>% 
  summarise(Resist=sum(Resistant),Sensi=sum(Sensitive),
            Susp=sum(Suspect),Tot=sum(Total))
dataCamem<-drop.levels(dataCamem)
data2map<-merge(dataCamem,coorddep,by.x="dptmt",by.y="dep_ID")

#mapping the results for each "programme"
colovec<-c(brewer.pal(9,"Reds")[7],brewer.pal(9,"Blues")[7],
           brewer.pal(9,"Oranges")[5])
for (i in 1:length(levels(data2map$themat_ID))){
  temp<-data2map[data2map$themat_ID==levels(data2map$themat_ID)[i],]
  png(file=paste("output/",temp$themat_ID,temp$pest,".png",sep=""),
      width=4,height=4,units="in",res=300)
  op<-par(mar=c(0,0,0,0))
  plot(departeLight,border="grey70")
  plot(regionsLight,lwd=2,add=TRUE)
  draw.pie(x=temp$longitude,y=temp$latitude,
           z=cbind((as.numeric(as.character(temp$Resist))),
                   (as.numeric(as.character(temp$Sensi))),
                   (as.numeric(as.character(temp$Susp)))),
           col=colovec,lty=0,
           radius=(sqrt(as.numeric(as.character(temp$Tot)))*15000),
           labels=NA)
  text(x=temp$longitude,y=temp$latitude,
       labels=as.character(temp$Tot),cex=1.2)
  par(op)
  dev.off()
}


##############################################################################/
#END
##############################################################################/
