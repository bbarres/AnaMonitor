##############################################################################/
##############################################################################/
#Script for Plasmopara viticola CI50 cyazo vs amisulbrom
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(drc)
library(plotrix)
library(gdata)
library(tidyr)
library(RColorBrewer)

#loading the data
datamyc<-read.table("data/2020-22_mildiouVigne_CDR_AmiCyaz.txt",
                     header=TRUE,stringsAsFactors=TRUE,sep=";")


##############################################################################/
#Regression analysis of mycelial growth experiment 14 days####
##############################################################################/

datamyc<-datamyc[datamyc$lect_echec!=1,]
datamyc<-drop.levels(datamyc)

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
CompRez<-data.frame(Subs_Act=factor(),sample_ID=factor(),read_time=factor(),
                    ED50=character(),ED95=character(),ED99=character())

#we make a subselection of the data according to the SA
pdf(file="output/plot_cyazo_vs_amisul.pdf",width=7)
for (j in 1:length(SAlist)) {
  data_subSA<-datamyc[datamyc$pest_sa_id==SAlist[j],]
  data_subSA$ech_id<-drop.levels(data_subSA$ech_id)
  
  REZSA<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                    read_time=factor(),ED50=character(),
                    ED95=character(),ED99=character())
  
  for (i in 1:dim(table(data_subSA$ech_id))[1]) {
    tempdat<-data_subSA[data_subSA$ech_id==names(table(data_subSA$ech_id))[i],]
    # if(tempdat[tempdat$dose==max(tempdat$dose),"rslt_03"]>50) {
    #   tempx<-data.frame("Subs_Act"=SAlist[j],"sample_ID"=tempdat$ech_id[1],
    #                     "read_time"=data_subSA$tps_expo[1],
    #                     "ED50"=paste(">",max(tempdat$dose),sep=""),
    #                     "ED95"=paste(">",max(tempdat$dose),sep=""),
    #                     "ED99"=paste(">",max(tempdat$dose),sep=""))
    # } else {
      temp.m1<-drm(rslt_03~dose,
                   data=tempdat,
                   fct=LL.3())
      plot(temp.m1,ylim=c(0,110),xlim=c(0,100),
           main=paste(SAlist[j],as.character(tempdat$ech_id[1])))
      temp<-ED(temp.m1,c(50,5,1),type="absolute")
      tempx<-data.frame("Subs_Act"=SAlist[j],
                        "sample_ID"=as.character(tempdat$ech_id[1]),
                        "read_time"=data_subSA$tps_expo[1],
                        "ED50"=as.character(temp[1]),
                        "ED95"=as.character(temp[2]),
                        "ED99"=as.character(temp[3]))
    # }
    
    REZSA<-rbind(REZSA,tempx)
  }
  CompRez<-rbind(CompRez,REZSA)
}
dev.off()

write.table(CompRez, file="output/results_cyazo_vs_amisul.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#Looking at the correlation between amisulbrom and cyazofamid CI50####
##############################################################################/

plot(log(as.numeric(CompRez$ED50[1:4])),log(as.numeric(CompRez$ED50[5:8])),
     xlab="CI50 amisulbrom SHAM",ylab="CI50 cyazofamid SHAM")

CompRez_wide<-spread(CompRez[,1:4],Subs_Act,ED50)
pdf(file="output/plot_correlAmisCyazo.pdf",width=7)
plot(as.numeric(CompRez_wide$AMISULBROM_SHAM),
     as.numeric(CompRez_wide$`CYAZOFAMIDE _SHAM`),
     xlab="CI50 amisulbrom SHAM",ylab="CI50 cyazofamid SHAM",
     xlog=TRUE,ylog=TRUE,log="xy",las=1)
dev.off()


##############################################################################/
#END
##############################################################################/