##############################################################################/
##############################################################################/
#Script for Plasmopara viticola 2023 CI50 oxathia and zoxa
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(drc)
library(plotrix)
library(gdata)
library(tidyr)
library(RColorBrewer)

#loading the data
datamyc<-read.table("data/2023_PdS_mildiouVigne_CDR_OxaZoxa.txt",
                    header=TRUE,stringsAsFactors=TRUE,sep=";")


##############################################################################/
#Regression analysis of sporulation  experiment ####
##############################################################################/

datamyc<-datamyc[datamyc$lect_echec!=1,]
datamyc<-drop.levels(datamyc)

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
CompRez<-data.frame(Subs_Act=factor(),sample_ID=factor(),read_time=factor(),
                    ED50=character(),ED95=character(),ED99=character())

#we make a subselection of the data according to the SA
pdf(file="output/2023_plot_DoseResp.pdf",width=7)
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

write.table(CompRez, file="output/2023_results_DosResp.txt",
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










##############################################################################/
##############################################################################/
#Analysis for Venturia sp by germination bioassay
##############################################################################/
##############################################################################/

source("creative_load_data.R")

#temporary data set loading waiting for the true final data set
VentuGer.dat<-creadat[creadat$species!="Alternaria sp" & 
                        creadat$species!="Alternaria alternata" &
                        creadat$species!="Alternaria arborescens" &
                        creadat$test_type=="spore_germ",]

#removing missing data
VentuGer.dat<-VentuGer.dat[!is.na(VentuGer.dat$perc_croiss),]
VentuGer.dat<-drop.levels(VentuGer.dat,reorder=FALSE)
#VentuGer.dat<-VentuGer.dat[!is.na(VentuGer.dat$perc_croiss),]

collist<-c("forestgreen","black")


##############################################################################/
#Analysis for Venturia sp by germination bioassay
##############################################################################/

REZ_VentuGe<-data.frame("strain_ID"=as.character(),"ActiveSub"=as.character(),
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

for (j in 1: length(levels(VentuGer.dat$active_substance))) {
  
  tempAS<-levels(VentuGer.dat$active_substance)[j]
  subdat<-VentuGer.dat[VentuGer.dat$active_substance==tempAS,]
  #some individual never reach an inhibition of 50%, event for the highest 
  #tested concentration. 
  subdat_rez<-as.character(subdat[subdat$dose==max(subdat$dose) & 
                                    subdat$perc_croiss>50,
                                  "strain_ID"])
  subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
  subdat_rez<-names(table(subdat_rez))
  ifelse(length(subdat_rez)==0,
         REZsub<-data.frame("strain_ID"=as.character(),
                            "ActiveSub"=as.character(),
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
                            "e-pval"=as.numeric()),
         REZsub<-data.frame("strain_ID"=subdat_rez,"ActiveSub"=tempAS,
                            "method"=levels(subdat$test_type),
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
  subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]
  subdat<-drop.levels(subdat,reorder=FALSE)
  #the subsequent analyses are unnecessary if all the strains are >maxdose
  if(dim(subdat)[1]==0) {
    REZ_VentuGe<-rbind(REZ_VentuGe,REZsub)
  } else {
    pdf(file=paste("output/VenGe_",tempAS,".pdf",sep=""),
        width=12,height=30)
    op<-par(mfrow=c(10,4))
    for (i in 1: dim(table(subdat$strain_ID))[1]) {
      datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
      couleur<-collist[as.numeric(datatemp$strain_type)]
      typeline<-as.numeric(datatemp$species)[1]
      print(as.character(datatemp$strain_ID[1]))
      if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
        tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                          "ActiveSub"=tempAS,
                          "method"=levels(subdat$test_type),
                          "ED50-abs"=NA,
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
        REZsub<-rbind(REZsub,tempx)
        plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
             main=names(table(subdat$strain_ID))[i])
      } else { tryCatch({
        temp.m1<-drm(perc_croiss~dose,
                     data=datatemp,
                     fct=LL.4())
        plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
             lty=typeline,type="confidence",
             main=names(table(subdat$strain_ID))[i])
        plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col="red",
             pch=4,type="obs",add=TRUE)
        plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
             lty=typeline,pch=19,add=TRUE)
      },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        if (!exists("temp.m1")){
          tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                            "ActiveSub"=tempAS,
                            "method"=levels(subdat$test_type),
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
          tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                            "ActiveSub"=tempAS,
                            "method"=levels(subdat$test_type),
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
        REZsub<-rbind(REZsub,tempx)
        rm(temp.m1)
      }
      
    }
    par(op)
    dev.off()
    
    REZ_VentuGe<-rbind(REZ_VentuGe,REZsub)
  }
  
}

write.table(REZ_VentuGe,file="output/Venturia_spore.txt",quote=FALSE,
            col.names=TRUE,row.names=FALSE,sep="\t")


##############################################################################/
#END
##############################################################################/




















##############################################################################/
##############################################################################/
#Analysis for Venturia by mycelial growth bioassay
##############################################################################/
##############################################################################/

source("creative_load_data.R")

#subsetting the global dataset
ventuGer.dat<-creadat[creadat$species=="V. inaequalis" | 
                        creadat$species=="V. asperata",]
ventuGer.dat<-ventuGer.dat[ventuGer.dat$test_type=="germination",]
ventuGer.dat$species<-factor(ventuGer.dat$species,
                             levels=rev(levels(ventuGer.dat$species)))
ventuGer.dat<-drop.levels(ventuGer.dat,reorder=FALSE)
#ventuGer.dat<-ventuGer.dat[!is.na(ventuGerc.dat$perc_croiss),]

collist<-c("forestgreen","black")


####boscalild
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="boscalid",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGeboscalid.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeBosc<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="boscalid")



####captane
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="captane",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGecaptane.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeCapt<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="captane")



####cyprodinil
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="cyprodinil",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGecyprodinil.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeCypr<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="cyprodinil")



####difenoconazole
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="difenoconazole",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGedifenoconazole.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeDife<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="difenoconazole")



####dithianon
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="dithianon",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGedithianon.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeDith<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="dithianon")



####dodine
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="dodine",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGedodine.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeDodi<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="dodine")



####mancozebe
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="mancozebe",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGemancozebe.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeManc<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="mancozebe")



####manebe
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="manebe",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGemanebe.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeMane<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="manebe")



####thirame
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="thirame",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGethirame.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeThir<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="thirame")



####trifloxystrobine
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="trifloxystrobine",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGetrifloxystrobine.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeTrif<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="trifloxystrobine")


#concatenating the different result files
REZ_VeGe<-rbind(REZ_VeGeBosc,REZ_VeGeCapt,REZ_VeGeCypr,REZ_VeGeDife,
                REZ_VeGeDith,REZ_VeGeDodi,REZ_VeGeManc,REZ_VeGeMane,
                REZ_VeGeThir,REZ_VeGeTrif)
write.table(REZ_VeGe,file="output/Venturia_germination.txt",quote=FALSE,
            col.names=TRUE,row.names=FALSE,sep="\t")


##############################################################################/
#END
##############################################################################/