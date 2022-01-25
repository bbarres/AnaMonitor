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
#one sample was not scored by two raters, so we remove it
compManipu<-compManipu[compManipu$ech_id!="20-044-02",]


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
write.table(CompRez,file="output/results_compManip.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#Statistical tests for comparison between raters####
##############################################################################/

#creating the empty result output file
CompRez2<-data.frame(Species=character(),Subs_Act=factor(),
                     sample_ID=factor(),
                     ED50.fr=character(),StErr.fr=character(),
                     ED50.ip=character(),StErr.ip=character(),
                     pval=character()
                     )

#we make a subselection of the data according to the SA
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
                           ED50.fr=character(),StErr.fr=character(),
                           ED50.ip=character(),StErr.ip=character(),
                           pval=character()
         ),
         REZSA<-data.frame("Species"=Esp_rez,
                           "Subs_Act"=SAlist[j],"sample_ID"=SA_rez,
                           "ED50.fr"=paste(">",max(data_subSA$dose),sep=""),
                           "StErr.fr"="no_eval",
                           "ED50.ip"=paste(">",max(data_subSA$dose),sep=""),
                           "StErr.ip"="no_eval",
                           "pval"="no_eval"
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
      temptest<-compParm(temp.m1,"e")
      temp<-ED(temp.m1,c(50),type="absolute")
      tempx<-data.frame("Species"=tempdat$species[1],
                        "Subs_Act"=SAlist[j],
                        "sample_ID"=tempdat$ech_id[1],
                        "ED50.fr"=as.character(temp[1,1]),
                        "StErr.fr"=as.character(temp[1,2]),
                        "ED50.ip"=as.character(temp[2,1]),
                        "StErr.ip"=as.character(temp[2,2]),
                        "pval"=as.character(temptest[4])
      )
      REZSA<-rbind(REZSA,tempx)}} else {
        REZSA<-REZSA
      }
  CompRez2<-rbind(CompRez2,REZSA)
}

#exporting the result as a text file
CompRez2<-CompRez2[order(CompRez2$Subs_Act,CompRez2$sample_ID),]
write.table(CompRez2,file="output/results_compManip2.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#END
##############################################################################/