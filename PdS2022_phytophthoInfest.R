##############################################################################/
##############################################################################/
#Analysis of Phytophthora infestans isolate of the 2022 monitoring
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(drc)
library(plotrix)
library(gdata)

phytophthoInfes<-read.table("data/2022_PdS_Phytophtora-infestans.txt",
                        header=TRUE,stringsAsFactors=TRUE,sep=";")


##############################################################################/
#Regression analysis of population germination for 2021 monitoring plan####
##############################################################################/

datamyc<-phytophthoInfes[phytophthoInfes$lect_echec!=1,]
datamyc<-drop.levels(datamyc)

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
#creating the empty result output file
CompRez<-data.frame("strain_ID"=as.character(),"ActiveSub"=as.character(),
                    "method"=as.character(),"ED50-abs"=as.numeric(),
                    "ED50-SE"=as.numeric(),"ED50-lower"=as.numeric(),
                    "ED50-upper"=as.numeric(),"b-param"=as.numeric(),
                    "b-SE"=as.numeric(),"b-tval"=as.numeric(),
                    "b-pval"=as.numeric(),"c-param"=as.numeric(),
                    "c-SE"=as.numeric(),"c-tval"=as.numeric(),
                    "c-pval"=as.numeric(),"d-param"=as.numeric(),
                    "d-SE"=as.numeric(),"d-tval"=as.numeric(),
                    "d-pval"=as.numeric(),"e-param"=as.numeric(),
                    "e-SE"=as.numeric(),"e-tval"=as.numeric(),
                    "e-pval"=as.numeric())

#we make a subselection of the data according to the SA
pdf(file="output/plot_phytophthoraInfes22.pdf",width=7)
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
         REZSA<-data.frame("strain_ID"=as.character(),
                           "ActiveSub"=as.character(),
                           "ED50-abs"=as.numeric(),
                           "ED50-SE"=as.numeric(),"ED50-lower"=as.numeric(),
                           "ED50-upper"=as.numeric(),"b-param"=as.numeric(),
                           "b-SE"=as.numeric(),"b-tval"=as.numeric(),
                           "b-pval"=as.numeric(),"c-param"=as.numeric(),
                           "c-SE"=as.numeric(),"c-tval"=as.numeric(),
                           "c-pval"=as.numeric(),"d-param"=as.numeric(),
                           "d-SE"=as.numeric(),"d-tval"=as.numeric(),
                           "d-pval"=as.numeric(),"e-param"=as.numeric(),
                           "e-SE"=as.numeric(),"e-tval"=as.numeric(),
                           "e-pval"=as.numeric()),
         REZSA<-data.frame("strain_ID"=subdat_rez,"ActiveSub"=SAlist[j],
                           "ED50-abs"=paste(">",max(subdat$dose),sep=""),
                           "ED50-SE"=NA,"ED50-lower"=NA,
                           "ED50-upper"=NA,"b-param"=NA,
                           "b-SE"=NA,"b-tval"=NA,
                           "b-pval"=NA,"c-param"=NA,
                           "c-SE"=NA,"c-tval"=NA,
                           "c-pval"=NA,"d-param"=NA,
                           "d-SE"=NA,"d-tval"=NA,
                           "d-pval"=NA,"e-param"=NA,
                           "e-SE"=NA,"e-tval"=NA,
                           "e-pval"=NA))
  
  #we limit the dataset to the sample that reach somehow a IC of 50%
  if(dim(data_subSA[!(data_subSA$ech_id %in% SA_rez),])[1]!=0) {
    SA.dat<-data_subSA[!(data_subSA$ech_id %in% SA_rez),]
    SA.dat<-drop.levels(SA.dat)
    for (i in 1:dim(table(SA.dat$ech_id))[1]) {
      tempdat<-SA.dat[SA.dat$ech_id==names(table(SA.dat$ech_id))[i],]
      tryCatch({
        temp.m1<-drm(rslt_03~dose,
                     data=tempdat,
                     fct=LL.4())
        plot(temp.m1,ylim=c(0,110),xlim=c(0,50),
             main=paste(tempdat$bioagr_id[1],
                        SAlist[j],names(table(SA.dat$ech_id))[i]))
      },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      if (!exists("temp.m1")){
        tempx<-data.frame("strain_ID"=names(table(data_subSA$ech_id))[i],
                          "ActiveSub"=SAlist[j],
                          "ED50-abs"="ERROR",
                          "ED50-SE"=NA,"ED50-lower"=NA,
                          "ED50-upper"=NA,"b-param"=NA,
                          "b-SE"=NA,"b-tval"=NA,
                          "b-pval"=NA,"c-param"=NA,
                          "c-SE"=NA,"c-tval"=NA,
                          "c-pval"=NA,"d-param"=NA,
                          "d-SE"=NA,"d-tval"=NA,
                          "d-pval"=NA,"e-param"=NA,
                          "e-SE"=NA,"e-tval"=NA,
                          "e-pval"=NA)
      } else {
        temp<-ED(temp.m1,50,type="absolute",interval="delta")
        tempb<-summary(temp.m1)$coefficients[1,]
        tempc<-summary(temp.m1)$coefficients[2,]
        tempd<-summary(temp.m1)$coefficients[3,]
        tempe<-summary(temp.m1)$coefficients[4,]
        tempx<-data.frame("strain_ID"=names(table(data_subSA$ech_id))[i],
                          "ActiveSub"=SAlist[j],
                          "ED50-abs"=temp[1],
                          "ED50-SE"=temp[2],"ED50-lower"=temp[3],
                          "ED50-upper"=temp[4],"b-param"=tempb[1],
                          "b-SE"=tempb[2],"b-tval"=tempb[3],
                          "b-pval"=tempb[4],"c-param"=tempc[1],
                          "c-SE"=tempc[2],"c-tval"=tempc[3],
                          "c-pval"=tempc[4],"d-param"=tempd[1],
                          "d-SE"=tempd[2],"d-tval"=tempd[3],
                          "d-pval"=tempd[4],"e-param"=tempe[1],
                          "e-SE"=tempe[2],"e-tval"=tempe[3],
                          "e-pval"=tempe[4])
      }
      REZSA<-rbind(REZSA,tempx)}} else {
        REZSA<-REZSA
      }
  CompRez<-rbind(CompRez,REZSA)
}
dev.off()

#exporting the result as a text file
CompRez<-CompRez[order(CompRez$ActiveSub,CompRez$strain_ID),]
write.table(CompRez, file="output/results_phytophthoraInfes22.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#END
##############################################################################/