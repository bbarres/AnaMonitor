##############################################################################/
##############################################################################/
#Analysis of Myzus spirotetramat bioassays
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(drc)
library(plotrix)
library(gdata)

dataReg<-read.table("data/ExtractionR_MyzusSpiro_72H.txt",
                    header=TRUE,stringsAsFactors=TRUE,sep=";")


##############################################################################/
#Observation of the raw data####
##############################################################################/

#aggregating the number of individuals per rep and dose
dataRegRep<-as.data.frame(aggregate(cbind(nb_vi,nb_mtot)~dose+dat_test,
                    data=dataReg,"sum"))
#the maximum number of individual for one repetition and one dose
ymasc<-max((dataRegRep$nb_vi+dataRegRep$nb_mtot))
#plotting the distribution of dead / alive for each repetition
op<-par(mfrow=c(2,2))
for (i in c(1:length(levels(as.factor(dataRegRep$dat_test))))) {
  temp<-dataRegRep[dataRegRep$dat_test==levels(as.factor(dataRegRep$dat_test))[i],]
  barplot(t(as.matrix(temp[,c(3:4)])),names.arg=temp$dose,
          las=1,xlab="Dose",ylab="Number",ylim=c(0,ymasc),
          main=levels(dataRegRep$dat_test)[i])
}
par(op)

#aggregating the number of individuals per dose
dataRegDos<-as.data.frame(aggregate(cbind(nb_vi,nb_mtot)~dose,
                                    data=dataReg,"sum"))
#barplot for all repetitions combined
barplot(t(as.matrix(dataRegDos[,c(2:3)])),
        names.arg=dataRegDos$dose,
        las=1,xlab="Dose",ylab="Number",
        main="All repetitions")


##############################################################################/
#Regression analysis of the first bioassay####
##############################################################################/

temp.m1<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
             weights=(nb_mtot+nb_vi),
             data=dataReg,
             curveid=dat_test,
             fct=LN.3u(),type="binomial")
plot(temp.m1,ylim=c(0,1.1),xlim=c(0,100),
     main=names(table(dataReg$ech_id))[1],
     ylab="mortality rate",legendPos=c(10,0.2))
compParm(temp.m1,"e")
EDcomp(temp.m1,c(50,50))
ED50m1<-ED(temp.m1,0.5,type="absolute")

temp.m2<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
             weights=(nb_mtot+nb_vi),
             data=dataRegRep,
             fct=LN.3u(),type="binomial")
ED50v<-ED(temp.m2,0.5,type="absolute")
plot(temp.m2,ylim=c(0,1.1),xlim=c(0,100),
     main=names(table(dataReg$ech_id))[1],
     type="confidence",ylab="mortality rate",
     las=1)
plot(temp.m2,ylim=c(0,1.1),xlim=c(0,100),
     main=names(table(dataReg$ech_id))[1],
     type="all",add=TRUE)
points((nb_mtot/(nb_mtot+nb_vi))~dose,data=dataRegDos,
       col="red",pch=19)
text(x=20,y=0.3,labels=paste("ED50=",round(ED50v[1],digits=4)))


##############################################################################/
#END
##############################################################################/