##############################################################################/
##############################################################################/
#Analysis of Venturia mycelial growth experiment on the 2020 monitoring
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(drc)
library(plotrix)
library(gdata)

#loading the data of the plan de surveillance
ventuCroi20<-read.table("data/2020_PdS_ventu_DOD.txt",header=TRUE,
                        stringsAsFactors=TRUE,sep=";")
#loading the data of the strain spores on film
SporeFilm20<-read.table("data/2020_results_ventuGerm_Film20.txt",header=TRUE,
                        stringsAsFactors=TRUE,sep=";")
SporeFilm20_2<-read.table("data/2020_results_ventuGerm_Film20_2.txt",header=TRUE,
                        stringsAsFactors=TRUE,sep=";")


##############################################################################/
#Regression analysis of population germination for 2020 monitoring plan####
##############################################################################/

datamyc<-ventuCroi20[ventuCroi20$lect_echec!=1,]

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
#creating the empty result output file
CompRez<-data.frame(Species=character(),Subs_Act=factor(),
                    sample_ID=factor(),read_time=factor(),
                    ED50=character(),ED95=character(),ED99=character())

#we make a subselection of the data according to the SA
pdf(file="output/plot_ventu_grow20.pdf",width=7)
for (j in 1:length(SAlist)) {
  data_subSA<-datamyc[datamyc$pest_sa_id==SAlist[j],]
  data_subSA$ech_id<-drop.levels(data_subSA$ech_id)
  #some individual never reach an inhibition of 50%, event for the highest 
  #tested concentration. 
  SA_rez<-as.character(data_subSA[data_subSA$dose==max(data_subSA$dose) 
                                  & data_subSA$rslt_03>50,
                                  "ech_id"])
  Esp_rez<-as.character(data_subSA[data_subSA$dose==max(data_subSA$dose) 
                                   & data_subSA$rslt_03>50,
                                   "bioagr_id"])
  ifelse(length(SA_rez)==0,
         REZSA<-data.frame(Species=character(),Subs_Act=factor(),
                           sample_ID=factor(),
                           read_time=factor(),ED50=character(),
                           ED95=character(),ED99=character()),
         REZSA<-data.frame("Species"=Esp_rez,
                           "Subs_Act"=SAlist[j],"sample_ID"=SA_rez,
                           "read_time"=data_subSA$tps_expo[1],
                           "ED50"=paste(">",max(data_subSA$dose),sep=""),
                           "ED95"=paste(">",max(data_subSA$dose),sep=""),
                           "ED99"=paste(">",max(data_subSA$dose),sep=""))
  )
  #we limit the data set to the sample that reach somehow a IC of 50%
  if(dim(data_subSA[!(data_subSA$ech_id %in% SA_rez),])[1]!=0) {
    SA.dat<-data_subSA[!(data_subSA$ech_id %in% SA_rez),]
    SA.dat<-drop.levels(SA.dat)
    for (i in 1:dim(table(SA.dat$ech_id))[1]) {
      tempdat<-SA.dat[SA.dat$ech_id==names(table(SA.dat$ech_id))[i],]
      temp.m1<-drm(rslt_03~dose,
                   data=tempdat,
                   fct=LL.4())
      plot(temp.m1,ylim=c(0,110),xlim=c(0,50),
           main=paste(data_subSA$bioagr_id[1],
                      SAlist[j],names(table(SA.dat$ech_id))[i]))
      temp<-ED(temp.m1,c(50,5,1),type="absolute")
      tempx<-data.frame("Species"=data_subSA$bioagr_id[1],
                        "Subs_Act"=SAlist[j],
                        "sample_ID"=names(table(SA.dat$ech_id))[i],
                        "read_time"=data_subSA$tps_expo[1],
                        "ED50"=as.character(temp[1]),
                        "ED95"=as.character(temp[2]),
                        "ED99"=as.character(temp[3]))
      REZSA<-rbind(REZSA,tempx)}} else {
        REZSA<-REZSA
      }
  CompRez<-rbind(CompRez,REZSA)
}
dev.off()

#exporting the result as a text file
CompRez<-CompRez[order(CompRez$Subs_Act,CompRez$sample_ID),]
write.table(CompRez, file="output/results_ventugrow20.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#Regression analysis of population germination for 2020 monitoring plan####
##############################################################################/

datamyc<-SporeFilm20[SporeFilm20_2$lect_echec!=1,]

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
#creating the empty result output file
CompRezFilm<-data.frame(Species=character(),Subs_Act=factor(),
                    sample_ID=factor(),read_time=factor(),
                    ED50=character(),ED95=character(),ED99=character())

#we make a subselection of the data according to the SA
pdf(file="output/plot_ventu_sporeFilm20.pdf",width=7)
for (j in 1:length(SAlist)) {
  data_subSA<-datamyc[datamyc$pest_sa_id==SAlist[j],]
  data_subSA$ech_id<-drop.levels(data_subSA$ech_id)
  #some individual never reach an inhibition of 50%, event for the highest 
  #tested concentration. 
  SA_rez<-as.character(data_subSA[data_subSA$dose==max(data_subSA$dose) 
                                  & data_subSA$rslt_03>50,
                                  "ech_id"])
  Esp_rez<-as.character(data_subSA[data_subSA$dose==max(data_subSA$dose) 
                                   & data_subSA$rslt_03>50,
                                   "bioagr_id"])
  ifelse(length(SA_rez)==0,
         REZSA<-data.frame(Species=character(),Subs_Act=factor(),
                           sample_ID=factor(),
                           read_time=factor(),ED50=character(),
                           ED95=character(),ED99=character()),
         REZSA<-data.frame("Species"=Esp_rez,
                           "Subs_Act"=SAlist[j],"sample_ID"=SA_rez,
                           "read_time"=data_subSA$tps_expo[1],
                           "ED50"=paste(">",max(data_subSA$dose),sep=""),
                           "ED95"=paste(">",max(data_subSA$dose),sep=""),
                           "ED99"=paste(">",max(data_subSA$dose),sep=""))
  )
  #we limit the data set to the sample that reach somehow a IC of 50%
  if(dim(data_subSA[!(data_subSA$ech_id %in% SA_rez),])[1]!=0) {
    SA.dat<-data_subSA[!(data_subSA$ech_id %in% SA_rez),]
    SA.dat<-drop.levels(SA.dat)
    for (i in 1:dim(table(SA.dat$ech_id))[1]) {
      tempdat<-SA.dat[SA.dat$ech_id==names(table(SA.dat$ech_id))[i],]
      temp.m1<-drm(rslt_03~dose,
                   data=tempdat,
                   fct=LL.4())
      plot(temp.m1,ylim=c(0,110),xlim=c(0,50),
           main=paste(data_subSA$bioagr_id[1],
                      SAlist[j],names(table(SA.dat$ech_id))[i]))
      temp<-ED(temp.m1,c(50,5,1),type="absolute")
      tempx<-data.frame("Species"=data_subSA$bioagr_id[1],
                        "Subs_Act"=SAlist[j],
                        "sample_ID"=names(table(SA.dat$ech_id))[i],
                        "read_time"=data_subSA$tps_expo[1],
                        "ED50"=as.character(temp[1]),
                        "ED95"=as.character(temp[2]),
                        "ED99"=as.character(temp[3]))
      REZSA<-rbind(REZSA,tempx)}} else {
        REZSA<-REZSA
      }
  CompRezFilm<-rbind(CompRezFilm,REZSA)
}
dev.off()

#exporting the result as a text file
CompRezFilm<-CompRezFilm[order(CompRezFilm$Subs_Act,CompRezFilm$sample_ID),]
write.table(CompRezFilm, file="output/results_ventuSporeFilm20.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#END
##############################################################################/