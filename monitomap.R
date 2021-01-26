##############################################################################/
##############################################################################/
#Different code for producing monitoring-related maps
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

#load geographical data of departements and regions
load("data/departeLight.RData")
load("data/regionsLight.RData")

#load the barycentre coordinates of departements and regions
load("data/coorddep.RData")
load("data/coordreg.RData")

#loading the data exported from PROSPER "tableau de bord" 2018
datasampl<-read.table("data/2018_toutheme_cor.txt",header=TRUE,
                      sep="\t",stringsAsFactors=TRUE)
datasampl<-merge(datasampl,coordreg,by.x="Region",by.y="reg_CODENAME")

#loading the data exported from PROSPER "tableau de bord" 2019
datasampl<-read.table("data/2019_prelevements.txt",header=TRUE,
                      sep="\t",stringsAsFactors=TRUE)
datasampl<-merge(datasampl,coordreg,by.x="Region",by.y="reg_CODENAME")

#loading the data exported from PROSPER "tableau de bord" 2020
datasampl<-read.table("data/2020_prelevements_BM.txt",header=TRUE,
                      sep="\t",stringsAsFactors=TRUE)
datasampl<-merge(datasampl,coordreg,by.x="Region",by.y="reg_CODENAME")

#loading the data for the flonicamid / Dysaphis experiment
datasampl<-read.table("data/floni_dysa_prelev.txt",header=TRUE,
                      sep="\t",stringsAsFactors=TRUE)
datasampl<-merge(datasampl,coordreg,by.x="Region",by.y="reg_CODENAME")

#loading the result data by departement
dataresult<-read.table("data/2017_themefin.txt",header=TRUE,sep="\t",
                       colClasses="character")
dataresult<-read.table("data/2018_themefin.txt",header=TRUE,sep="\t",
                       colClasses="character")
dataresult<-read.table("data/2019_rezultBM.txt",header=TRUE,sep="\t",
                       colClasses="character")
dataresult<-read.table("data/floni_dysa_rez.txt",header=TRUE,sep="\t",
                       colClasses="character")
dataresult<-read.table("data/2020_mild_rez.txt",header=TRUE,sep="\t",
                       colClasses="character")


#to streamline subsequent analysis, we turned the resistance status factor
#into two different columns
dataresult$rslt_RS<-factor(dataresult$rslt_RS,levels=c("R","S"))
dataresult$Resistant<-as.numeric(dataresult$rslt_RS==
                                   levels(dataresult$rslt_RS)[1])
dataresult$Sensitive<-as.numeric(dataresult$rslt_RS==
                                   levels(dataresult$rslt_RS)[2])
dataresult$Total<-dataresult$Resistant+dataresult$Sensitive


##############################################################################/
#Monitoring sampling maps by regions####
##############################################################################/

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


##############################################################################/
#Maps of the results of the monitoring####
##############################################################################/

#a note for futur-me: this code will probably need some minor updates in order
#to work properly on the raw exported-from-PROSPER file

#grouping and counting the resistant and sensitive samples
dataresult$themat_ID<-paste(dataresult$themat_ID,dataresult$SA)
dataresult$themat_ID<-as.factor(dataresult$themat_ID)
dataCamem<-dataresult %>% 
  group_by(themat_ID,dptmt,pest,host) %>% 
  summarise(Resist=sum(Resistant),Sensi=sum(Sensitive),Tot=sum(Total))
data2map<-merge(dataCamem,coorddep,by.x="dptmt",by.y="dep_ID")

#mapping the results for each "programme" of the monitoring
colovec<-c(brewer.pal(9,"Reds")[7],brewer.pal(9,"Blues")[7])
for (i in 1:length(levels(data2map$themat_ID))){
  temp<-data2map[data2map$themat_ID==levels(data2map$themat_ID)[i],]
  png(file=paste("output/",temp$themat_ID,temp$pest,".png",sep=""),
      width=4,height=4,units="in",res=300)
  op<-par(mar=c(0,0,0,0))
  plot(departeLight,border="grey70")
  plot(regionsLight,lwd=2,add=TRUE)
  draw.pie(x=temp$longitude,y=temp$latitude,
           z=cbind((as.numeric(as.character(temp$Resist))),
                   (as.numeric(as.character(temp$Sensi)))),
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