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


#aggregating the number of individuals per rep and dose
dataRegRep<-as.data.frame(aggregate(cbind(nb_vi,nb_mtot)~dose+rep_test,
                    data=dataReg,"sum"))

#aggregating the number of individuals per dose
dataRegDos<-as.data.frame(aggregate(cbind(nb_vi,nb_mtot)~dose,
                                    data=dataReg,"sum"))

#the maximum number of individual for one repetition and one dose
ymasc<-max((dataRegRep$nb_vi+dataRegRep$nb_mtot))
#plotting the distribution of dead / alive for each repetition
op<-par(mfrow=c(3,4))
for (i in c(1:length(levels(dataRegRep$rep_test)))) {
  temp<-dataRegRep[dataRegRep$rep_test==levels(dataRegRep$rep_test)[i],]
  barplot(t(as.matrix(temp[,c(3:4)])),names.arg=temp$dose,
          las=1,xlab="Dose",ylab="Number",ylim=c(0,ymasc),
          main=levels(dataRegRep$rep_test)[i])
}
par(op)


##############################################################################/
#Regression analysis of the first bioassay####
##############################################################################/

temp.m1<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
             weights=(nb_mtot+nb_vi),
             data=dataReg,
             curveid=rep_test,
             fct=LN.3u(),type="binomial")
plot(temp.m1,ylim=c(0,1.1),xlim=c(0,100),
     main=names(table(dataReg$ech_id))[1],
     ylab="mortality rate")
compParm(temp.m1,"e")

temp.m2<-drm(nb_mtot/(nb_mtot+nb_vi)~dose,
             weights=(nb_mtot+nb_vi),
             data=dataReg,
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