#small script to extract consumption data from NHS - example October 2017

#CHANGE this to your working directory (where the files are stored)
setwd("C:/Users/Brett Sallach/Desktop/R Stuff")

#load the raw data
data_raw <- read.csv("1512_Detailed_Prescribing_Information.csv",header = TRUE, stringsAsFactors = FALSE, sep = ",")
practices <- read.csv("Practices.csv",header=FALSE,stringsAsFactors = FALSE, sep=",")
chems <- read.csv("API_BNF_Penicillins.csv",header=FALSE,stringsAsFactors = FALSE, sep=",")

#select only the relevant practices
data_York <- data_raw[data_raw$Practice.Name%in%practices$V1,]

#add column indicating the API code (first 9 characters of BNF code)
data_York$API <- substr(data_York$BNF.Code,start=1,stop=9)

#select only the relevant APIs (based on BNF codes)
data_York <- data_York[data_York$API%in%chems$V1,]

#create list per relevant API
data_York_sub <- vector(mode = "list", length = length(chems$V1))
for (i in 1:length(data_York_sub)) {
  data_York_sub[[i]] <- data_York[data_York$API==chems$V1[i],]
}

#Compute total dose per API (simple formulations where 1 value is mentioned)
for (i in 1:length(data_York_sub)) {
  if (nrow(data_York_sub[[i]])>0) {
    for (j in 1:nrow(data_York_sub[[i]])) {
      #check whether multiple doses are mentioned in the description
      data_York_sub[[i]]$n_doses[j] <- length(as.numeric(unlist(regmatches(data_York_sub[[i]]$BNF.Description[j],gregexpr("[0-9]+",data_York_sub[[i]]$BNF.Description[j])))))
      #add value for dosage only for those records where one specific dosage is mentioned
      data_York_sub[[i]]$dose[j] <- ifelse(data_York_sub[[i]]$n_items[j]==1,
                                           as.numeric(unlist(regmatches(data_York_sub[[i]]$BNF.Description[j],gregexpr("[0-9]+",data_York_sub[[i]]$BNF.Description[j])))),
                                           NA)
      #compute the total sum for these records
      data_York_sub[[i]]$total_dose[j] <- data_York_sub[[i]]$dose[j] * data_York_sub[[i]]$Quantity[j]
      
    }
    #write the relevant APIs (BNF codes) to your directory
    write.csv(data_York_sub[[i]],file=paste(Sys.Date(),"_",chems$V1[i],".csv",sep=""),row.names = FALSE)
    }
  } 