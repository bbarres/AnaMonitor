##############################################################################/
##############################################################################/
#Comparison of same experiment scored by two persons
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(drc)
library(plotrix)
library(gdata)

#loading the data
compManipu<-read.table("data/ventu_manipul.txt",header=TRUE,
                      stringsAsFactors=TRUE,sep="\t")


##############################################################################/
#Regression analysis of germination tube length####
##############################################################################/

datamyc<-compManipu
datamyc$rslt_03<-as.numeric(as.character(datamyc$rslt_03))

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
#creating the empty result output file
CompRez<-data.frame(Species=character(),Subs_Act=factor(),
                    sample_ID=factor(),repetition=factor(),
                    ED50=character(),StErr=character())

#we make a subselection of the data according to the SA
pdf(file="output/plot_compManip.pdf",width=7)
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
  #we limit the data set to the sample that reach somehow a IC of 50%
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
write.table(CompRez, file="output/results_compManip.txt",
            sep="\t",quote=FALSE,row.names=FALSE)

#creating the empty result output file
CompRez<-data.frame(Species=character(),Subs_Act=factor(),
                    sample_ID=factor(),repetition=factor(),
                    ED50=character(),StErr=character())

#we make a subselection of the data according to the SA
pdf(file="output/plot_compManip.pdf",width=7)
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
  #we limit the data set to the sample that reach somehow a IC of 50%
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
write.table(CompRez, file="output/results_compManip.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#END
##############################################################################/