##############################################################################/
##############################################################################/
#Script to produce the JER2024 figures
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
rezMoni<-read.table("data/PlasmoVitiJER2024.txt",header=TRUE,sep="\t",
                    colClasses="character")
#turning the resistance status factor into two different columns
rezMoni$rslt_RS<-factor(rezMoni$rslt_RS,levels=c("R","S"))
rezMoni$Resistant<-as.numeric(rezMoni$rslt_RS==levels(rezMoni$rslt_RS)[1])
rezMoni$Sensitive<-as.numeric(rezMoni$rslt_RS==levels(rezMoni$rslt_RS)[2])
rezMoni$Total<-rezMoni$Resistant+rezMoni$Sensitive
#grouping and counting the resistant and sensitive samples
rezMoni$themat_ID<-(rezMoni$SA)
rezMoni$themat_ID<-as.factor(rezMoni$themat_ID)
dataCamem<-rezMoni %>% 
  group_by(themat_ID,dptmt) %>% 
  summarise(Resist=sum(Resistant),Sensi=sum(Sensitive),Tot=sum(Total))
data2map<-merge(dataCamem,coorddep,by.x="dptmt",by.y="dep_ID")

#defining a set of colors
colovec<-c(brewer.pal(9,"Reds")[7],brewer.pal(9,"Blues")[7])


##############################################################################/
#Maps of the results of the monitoring cumulated between 2018-2023####
##############################################################################/

#mapping the results for each "programme"
colovec<-c(brewer.pal(9,"Reds")[7],brewer.pal(9,"Blues")[7])
for (i in 1:length(levels(data2map$themat_ID))){
  temp<-data2map[data2map$themat_ID==levels(data2map$themat_ID)[i],]
  png(file=paste("output/","JER_",temp$themat_ID,temp$pest,".png",sep=""),
      width=4,height=4,units="in",res=300)
  op<-par(mar=c(0,0,0,0))
  plot(DEP_SHP,border="grey70")
  plot(REG_SHP,lwd=2,add=TRUE)
  draw.pie(x=temp$longitude,y=temp$latitude,
           z=cbind((as.numeric(as.character(temp$Resist))),
                   (as.numeric(as.character(temp$Sensi)))),
           col=colovec,lty=0,
           radius=(sqrt(as.numeric(as.character(temp$Tot)))*9000),
           labels=NA)
  text(x=temp$longitude,y=temp$latitude,
       labels=as.character(temp$Tot),cex=1.2)
  #scalebar(c(191260,6060000),300000,"km",division.cex=0.8)
  par(op)
  dev.off()
}


##############################################################################/
#barplot of occurence evolution####
##############################################################################/

#grouping and counting the resistant and sensitive samples by year
dataYear<-rezMoni %>% 
  group_by(themat_ID,year) %>% 
  summarise(Resist=sum(Resistant),Sensi=sum(Sensitive),Tot=sum(Total))

colovec<-c(brewer.pal(9,"Reds")[7],brewer.pal(9,"Blues")[7])

#plotting
for (i in 1:length(levels(dataYear$themat_ID))){
temp1<-dataYear[dataYear$themat_ID==levels(dataYear$themat_ID)[i],]
nametemp<-levels(dataYear$themat_ID)[i]
yeartemp<-temp1$year
temp1<-prop.table(as.matrix(temp1[,3:4]),margin=1)*100
row.names(temp1)<-yeartemp

png(file=paste("output/","occurence_",nametemp,".png",sep=""),
    width=8,height=6,units="in",res=300)
barplot(t(temp1),col=colovec,las=1,main=nametemp,cex.names=1.5,
        font=2,space=0.7)
dev.off()}


##############################################################################/
#Distribution of frequencies within population####
##############################################################################/

#distribution between 2018-2023
for (i in 1:length(levels(rezMoni$themat_ID))){
  temp<-rezMoni[rezMoni$themat_ID==levels(rezMoni$themat_ID)[i],]
  png(file=paste("output/","2023_hist",temp$themat_ID,temp$pest,".png",sep=""),
      width=4,height=7,units="in",res=300)
  hist(as.numeric(temp$Perc_spor),ann=FALSE,las=1,freq=TRUE,
       include.lowest=TRUE,breaks=seq(0,100,by=10))
  dev.off()
}

#distribution for each year
#loading the data exported from PROSPER "tableau de bord"
rezMoni$themaYear<-paste(rezMoni$SA,rezMoni$year,sep="_")
rezMoni$themaYear<-as.factor(rezMoni$themaYear)

#plot histogram for the distribution of the percentage of sporulation by year
for (i in 1:length(levels(rezMoni$themaYear))){
  temp<-rezMoni[rezMoni$themaYear==levels(rezMoni$themaYear)[i],]
  png(file=paste("output/","JER_hist",temp$themaYear,".png",sep=""),
      width=5,height=4,units="in",res=300)
  hist(as.numeric(temp$Perc_spor),ann=FALSE,las=1,freq=TRUE,
       include.lowest=TRUE,breaks=seq(0,100,by=10),
       xlim=c(0,100))
  abline(v=mean(as.numeric(temp$Perc_spor)),lty=2,col="red",lwd=3)
  dev.off()
}


#evolution of mean frequencies in R population
temp<-rezMoni[rezMoni$rslt_RS=="R",]
temp<-rezMoni
temp<-droplevels(temp)
temp1<-aggregate(as.numeric(temp$Perc_spor),list(temp$themaYear),FUN=mean)
temp1$SA<-c(rep("ametoc",6),c(rep("amisu",6)),c(rep("AOX",6)),
            c(rep("cyazo",6)),c(rep("fluop",5)),c(rep("rien",8)))

write.table(temp1,file="output/JER2024meanFreq.txt",sep="\t",
            quote=FALSE)

plot(temp1[temp1$SA=="ametoc",]$x,type="b",las=1,
     ylab="Fréquence",xlab="Année")


##############################################################################/
#END
##############################################################################/