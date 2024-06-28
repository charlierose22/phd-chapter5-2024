#small script to extract consumption data from NHS - example October 2017
cat("\014"); rm(list=ls())
st=Sys.time(); setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load package xlsx to write to excel spreadsheet
# install.packages("xlsx") #install required package, only needed once
require("xlsx")

#load the raw data
sheets <- c(1512,1601:1609)
class <- "UTI"

for (p in sheets) {
flnm <- paste0(p,"_Detailed_Prescribing_Information.csv")

data_raw <- read.csv(flnm,header = TRUE, stringsAsFactors = FALSE, sep = ",")
practices <- read.csv("Practices.csv",header=FALSE,stringsAsFactors = FALSE, sep=",")
chems <- read.csv(paste0("API_BNF_",class,".csv"),header=FALSE,stringsAsFactors = FALSE, sep=",")

#select only the relevant practices
data_York <- data_raw[data_raw$Practice.Name%in%practices$V1,]

#add column indicating the API code (first 9 characters of BNF code)
data_York$API <- substr(data_York$BNF.Code,start=1,stop=9)

#select only the relevant APIs (based on BNF codes)
data_York <- data_York[data_York$API%in%chems$V1,]
data_York$BNF.Description[data_York$BNF.Description=="Augmentin_Pdr For Susp 250/62mg/5ml S/F"] <- "Augmentin_Pdr For Susp 250mg/62mg/5ml S/F"
data_York$BNF.Description[data_York$BNF.Description=="Augmentin_Pdr For Susp 125/31mg/5ml S/F"] <- "Augmentin_Pdr For Susp 125mg/62mg/5ml S/F"

data_York$BNF.Description[data_York$BNF.Description=="Co-Trimoxazole_Oral Susp 40/200mg/5mlS/F"] <- "Co-Trimoxazole_Oral Susp 40mg/200mg/5mlS/F"

data_York$BNF.Description[data_York$BNF.Description=="Tetralysal 300_Cap"] <- "Tetralysal 300mg_Cap"

#create list per relevant API
data_York_sub <- vector(mode = "list", length = length(chems$V1))

for (i in 1:length(data_York_sub)) {
  data_York_sub[[i]] <- data_York[data_York$API==chems$V1[i],]
  
  if (nrow(data_York_sub[[i]])>0) {
    for (j in 1:nrow(data_York_sub[[i]])) {
      data_York_sub[[i]]$n_doses[j] <- length(as.numeric(unlist(regmatches(data_York_sub[[i]]$BNF.Description[j],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",data_York_sub[[i]]$BNF.Description[j])))))
    }
    dosemax <- max(data_York_sub[[i]]$n_doses)
    
    #create empty columns for doses, units, and total
    data_York_sub[[i]][paste(rep("dose",dosemax),seq(1,dosemax))] <- NA
    data_York_sub[[i]][paste(rep("units",dosemax),seq(1,dosemax))] <- NA
    data_York_sub[[i]][paste(rep("total",dosemax),seq(1,dosemax))] <- NA
    
    #fill columns with doses and with units
    for (j in 1:nrow(data_York_sub[[i]])) {
      data_York_sub[[i]][j,paste(rep("dose",data_York_sub[[i]]$n_doses[j]),seq(1,data_York_sub[[i]]$n_doses[j]))] <-
        as.numeric(unlist(regmatches(data_York_sub[[i]]$BNF.Description[j],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",data_York_sub[[i]]$BNF.Description[j]))))
      
      data_York_sub[[i]][j,paste(rep("units",data_York_sub[[i]]$n_doses[j]),seq(1,data_York_sub[[i]]$n_doses[j]))] <-
        substr(unlist(strsplit(data_York_sub[[i]]$BNF.Description[j], "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE))[-1],start=1,stop=2)
      
      data_York_sub[[i]][j,paste(rep("total",data_York_sub[[i]]$n_doses[j]),seq(1,data_York_sub[[i]]$n_doses[j]))] <-
        as.numeric(unlist(regmatches(data_York_sub[[i]]$BNF.Description[j],gregexpr("[[:digit:]]+\\.*[[:digit:]]*",data_York_sub[[i]]$BNF.Description[j]))))*
        data_York_sub[[i]]$Quantity[j]*data_York_sub[[i]]$Items[j]
      
      
    }
    
    #write to worksheet  
    if (i == 1) {
      write.xlsx(x = data_York_sub[[i]], file = paste0(Sys.Date(),class,p,".xlsx"), sheetName = chems$V1[i], row.names = FALSE)
    } else {
      write.xlsx(x = data_York_sub[[i]], file = paste(Sys.Date(),class,p,".xlsx",sep=""), sheetName = chems$V1[i], append = TRUE, row.names = FALSE)
    }
  }
  

}
}
