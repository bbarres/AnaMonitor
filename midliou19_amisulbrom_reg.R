##############################################################################/
##############################################################################/
#Script for milidou 2019 amisulbrom results
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(drc)
library(plotrix)
library(gdata)
library(tidyr)
library(RColorBrewer)

#loading the data
datamyc2<-read.table("data/2019_VEGEmildiou_rez.txt",header=TRUE,
                     stringsAsFactors=TRUE,sep=";")


##############################################################################/
#Regression analysis of mycelial growth experiment 14 days####
##############################################################################/

datamyc<-datamyc2[datamyc2$pest_sa_id=="AMISULBROM_SHAM" & 
                    datamyc2$meth_id!="SPORUL_PLASMOPARA_VITICOLA_2018",]
#because it cause a crash of the loop, we remove one analyze --> to be fixed
datamyc<-datamyc[datamyc$ana_id!=list(1144,1077),]
datamyc<-drop.levels(datamyc)

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
CompRez<-data.frame(Subs_Act=factor(),sample_ID=factor(),read_time=factor(),
                    ED50=character(),ED95=character(),ED99=character())

#we make a subselection of the data according to the SA
pdf(file="output/plot_Amisulbrom19.pdf",width=7)
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
      plot(temp.m1,type="obs",add=TRUE,col="red")
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

write.table(CompRez, file="output/results_amislbrom19.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#END
##############################################################################/