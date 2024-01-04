##############################################################################/
##############################################################################/
#Script for pyrenophora in 2022 and 2023 CI50
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(drc)
library(plotrix)
library(gdata)
library(tidyr)
library(RColorBrewer)

#loading the data
datamyc<-read.table("data/2023_PdS_pyrenophora_2022-23.txt",
                    header=TRUE,stringsAsFactors=TRUE,sep=";")


##############################################################################/
#Regression analysis of mycelial growth  experiment####
##############################################################################/

datamyc<-datamyc[datamyc$lect_echec!=1,]
datamyc<-drop.levels(datamyc)

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
CompRez<-data.frame("strain_ID"=factor(),"ActiveSub"=factor(),
                    "ED50-abs"=as.numeric(),"ED50-SE"=as.numeric(),
                    "ED50-lower"=as.numeric(),"ED50-upper"=as.numeric(),
                    "b-param"=as.numeric(),"b-SE"=as.numeric(),
                    "b-tval"=as.numeric(),"b-pval"=as.numeric(),
                    "d-param"=as.numeric(),"d-SE"=as.numeric(),
                    "d-tval"=as.numeric(),"d-pval"=as.numeric(),
                    "e-param"=as.numeric(),"e-SE"=as.numeric(),
                    "e-tval"=as.numeric(),"e-pval"=as.numeric())

#we make a subselection of the data according to the SA
pdf(file="output/2023_PdSplot_pyreno.pdf",width=7)
for (j in 1:length(SAlist)) {
  data_subSA<-datamyc[datamyc$pest_sa_id==SAlist[j],]
  data_subSA$ech_id<-drop.levels(data_subSA$ech_id)
  
  REZSA<-data.frame("strain_ID"=factor(),"ActiveSub"=factor(),
                    "ED50-abs"=as.numeric(),"ED50-SE"=as.numeric(),
                    "ED50-lower"=as.numeric(),"ED50-upper"=as.numeric(),
                    "b-param"=as.numeric(),"b-SE"=as.numeric(),
                    "b-tval"=as.numeric(),"b-pval"=as.numeric(),
                    "d-param"=as.numeric(),"d-SE"=as.numeric(),
                    "d-tval"=as.numeric(),"d-pval"=as.numeric(),
                    "e-param"=as.numeric(),"e-SE"=as.numeric(),
                    "e-tval"=as.numeric(),"e-pval"=as.numeric())
  
  for (i in 1:dim(table(data_subSA$ech_id))[1]) {
    tempdat<-data_subSA[data_subSA$ech_id==names(table(data_subSA$ech_id))[i],]
    if(any(tempdat[tempdat$dose==max(tempdat$dose),"rslt_03"]>50)) {
      tempx<-data.frame("strain_ID"=tempdat$ech_id[1],"ActiveSub"=SAlist[j],
                        "ED50-abs"=paste(">",max(tempdat$dose),sep=""),
                        "ED50-SE"=NA,"ED50-lower"=NA,"ED50-upper"=NA,
                        "b-param"=NA,"b-SE"=NA,"b-tval"=NA,"b-pval"=NA,
                        "d-param"=NA,"d-SE"=NA,"d-tval"=NA,"d-pval"=NA,
                        "e-param"=NA,"e-SE"=NA,"e-tval"=NA,"e-pval"=NA)
    } else { tryCatch({
      temp.m1<-drm(rslt_03~dose,
                   data=tempdat,
                   fct=LL.3())
      plot(temp.m1,ylim=c(-10,120),xlim=c(0,100),type="confidence",
           main=paste(SAlist[j],as.character(tempdat$ech_id[1])))
      plot(temp.m1,ylim=c(-10,120),xlim=c(0,100),col="red",
           pch=4,type="obs",add=TRUE)
      
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      if (!exists("temp.m1")){
        tempx<-data.frame("strain_ID"=tempdat$ech_id[1],"ActiveSub"=SAlist[j],
                          "ED50-abs"="ERROR",
                          "ED50-SE"=NA,"ED50-lower"=NA,"ED50-upper"=NA,
                          "b-param"=NA,"b-SE"=NA,"b-tval"=NA,"b-pval"=NA,
                          "d-param"=NA,"d-SE"=NA,"d-tval"=NA,"d-pval"=NA,
                          "e-param"=NA,"e-SE"=NA,"e-tval"=NA,"e-pval"=NA)
      } else {
        temp<-ED(temp.m1,50,type="absolute",interval="delta")
        tempb<-summary(temp.m1)$coefficients[1,]
        tempd<-summary(temp.m1)$coefficients[2,]
        tempe<-summary(temp.m1)$coefficients[3,]
        tempx<-data.frame("strain_ID"=tempdat$ech_id[1],"ActiveSub"=SAlist[j],
                          "ED50-abs"=temp[1],"ED50-SE"=temp[2],
                          "ED50-lower"=temp[3],"ED50-upper"=temp[4],
                          "b-param"=tempb[1],"b-SE"=tempb[2],
                          "b-tval"=tempb[3],"b-pval"=tempb[4],
                          "d-param"=tempd[1],"d-SE"=tempd[2],
                          "d-tval"=tempd[3],"d-pval"=tempd[4],
                          "e-param"=tempe[1],"e-SE"=tempe[2],
                          "e-tval"=tempe[3],"e-pval"=tempe[4])
      }
    }
    REZSA<-rbind(REZSA,tempx)
    rm(temp.m1)
  }
  CompRez<-rbind(CompRez,REZSA)
}
dev.off()

write.table(CompRez,file="output/2023_PdSresults_pyreno.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#END
##############################################################################/