##############################################################################/
##############################################################################/
#Analysis of Alternaria mycelial growth experiment on the 2021 monitoring
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(drc)
library(plotrix)
library(gdata)

#loading the data
alterna21<-read.table("data/data_alternaria.txt",header=TRUE,
                        stringsAsFactors=TRUE,sep="\t")


##############################################################################/
#Regression analysis of mycelial growth for 2021 monitoring plan####
##############################################################################/

datamyc<-alterna21[alterna21$rslt_03!="echec",]
datamyc$rslt_03<-as.numeric(as.character(datamyc$rslt_03))

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
#creating the empty result output file
CompRez<-data.frame(Species=character(),Subs_Act=factor(),
                    sample_ID=factor(),repetition=factor(),
                    ED50=character())

#we make a subselection of the data according to the SA
pdf(file="output/plot_alter_21.pdf",width=7)
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
                                   "species"])
  ifelse(length(SA_rez)==0,
         REZSA<-data.frame(Species=character(),Subs_Act=factor(),
                           sample_ID=factor(),
                           repetition=factor(),ED50=character()
                           ),
         REZSA<-data.frame("Species"=Esp_rez,
                           "Subs_Act"=SAlist[j],"sample_ID"=SA_rez,
                           "repetition"="pooled",
                           "ED50"=paste(">",max(data_subSA$dose),sep="")
                           )
  )
  #we limit the dataset to the sample that reach somehow a IC of 50%
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
      plot(temp.m1,ylim=c(0,110),xlim=c(0,50),
           main=paste(data_subSA$bioagr_id[1],
                      SAlist[j],names(table(SA.dat$ech_id))[i]),
           type="obs",col="red",add=TRUE)
      temp<-ED(temp.m1,c(50),type="absolute")
      tempx<-data.frame("Species"=tempdat$species[1],
                        "Subs_Act"=SAlist[j],
                        "sample_ID"=tempdat$ech_id[1],
                        "repetition"="pooled",
                        "ED50"=as.character(temp[1])
                        )
      REZSA<-rbind(REZSA,tempx)}} else {
        REZSA<-REZSA
      }
  CompRez<-rbind(CompRez,REZSA)
}
dev.off()

#exporting the result as a text file
CompRez<-CompRez[order(CompRez$Subs_Act,CompRez$sample_ID),]
write.table(CompRez, file="output/results_alterna21.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#Regression analysis of mycelial growth by repet for 2021 monitoring plan####
##############################################################################/

datamyc<-alterna21[alterna21$rslt_03!="echec",]
datamyc$rslt_03<-as.numeric(as.character(datamyc$rslt_03))

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
#creating the empty result output file
CompRez<-data.frame(Species=character(),Subs_Act=factor(),
                    sample_ID=factor(),repetition=factor(),
                    ED50=character(),StErr=character())

#we make a subselection of the data according to the SA
pdf(file="output/plot_alter_21_repet.pdf",width=7)
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
                                   "species"])
  Rep_rez<-as.character(data_subSA[data_subSA$dose==max(data_subSA$dose) 
                                   & data_subSA$rslt_03>50,
                                   "repet"])
  ifelse(length(SA_rez)==0,
         REZSA<-data.frame(Species=character(),Subs_Act=factor(),
                           sample_ID=factor(),repetition=factor(),
                           ED50=character(),StErr=character()
         ),
         REZSA<-data.frame("Species"=Esp_rez,
                           "Subs_Act"=SAlist[j],"sample_ID"=SA_rez,
                           "repetition"=Rep_rez,
                           "ED50"=paste(">",max(data_subSA$dose),sep=""),
                           "StErr"="no_eval"
         )
  )
  #we limit the dataset to the sample that reach somehow a IC of 50%
  if(dim(data_subSA[!(data_subSA$ech_id %in% SA_rez),])[1]!=0) {
    SA.dat<-data_subSA[!(data_subSA$ech_id %in% SA_rez),]
    SA.dat<-drop.levels(SA.dat)
    for (i in 1:dim(table(SA.dat$ech_id))[1]) {
      tempdat<-SA.dat[SA.dat$ech_id==names(table(SA.dat$ech_id))[i],]
      temp.m1<-drm(rslt_03~dose,
                   data=tempdat,
                   curveid=repet,
                   fct=LL.4())
      plot(temp.m1,ylim=c(0,110),xlim=c(0,50),
           main=paste(data_subSA$bioagr_id[1],
                      SAlist[j],names(table(SA.dat$ech_id))[i]))
      temp<-ED(temp.m1,c(50),type="absolute")
      tempx<-data.frame("Species"=tempdat$species[1],
                        "Subs_Act"=SAlist[j],
                        "sample_ID"=tempdat$ech_id[1],
                        "repetition"=dimnames(temp)[[1]],
                        "ED50"=as.character(temp[,1]),
                        "StErr"=as.character(temp[,2])
      )
      REZSA<-rbind(REZSA,tempx)}} else {
        REZSA<-REZSA
      }
  CompRez<-rbind(CompRez,REZSA)
}
dev.off()

#exporting the result as a text file
CompRez<-CompRez[order(CompRez$Subs_Act,CompRez$sample_ID),]
write.table(CompRez, file="output/results_alterna21_repet.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#END
##############################################################################/