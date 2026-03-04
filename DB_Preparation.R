setwd("~/Dropbox/PONDERFUL/Database/data/Ecosystem functions/Analisis/All/SEM/Article/R1_11Sep2025/R1_documents/To Submit/DataScript")
source("~/Dropbox/PONDERFUL/Database/data/Ecosystem functions/Analisis/All/SEM/Article/Documents to submit/Database_Script/PackagesFuntions.R")

#set up your tables
EF_SEM <- read_excel("EMF_rawData.xlsx", sheet = "FullDB")

EF <- EF_SEM #N=206
PondID_Country <- EF[ , c(1,2)] #extract pond id and country
EF <- column_to_rownames(EF, var="Pond_ID")
EF <- EF[ , -c(1,20,21)] #extract the binomial variables

###****** Permutate the NA of your functions and taxonomic group 
EMF = EF[,c(1:5)]
EMF <- mice(EMF, method = "pmm", m = 5)  # Permutation with m = 5. Increase if necessary
EMF <- complete(EMF)  # extract the data
#Transform your EF to reach normatily
EMF$logCH4 <- log(EMF$CH4) 
EMF$logCHLa <- log(EMF$CHLa)
EMF$logBMzp <- log(EMF$BMzp)
EMF$logBMmi <- log(EMF$BMmi)
EMF = EMF[,c(6:9,5)]

###****** Calculate multidiversity
#This is calculated in the same way that weight EMF from Wang 2019
MDiv = EF[,c(6:9)]
max_richness <- apply(MDiv, 2, max, na.rm = TRUE)
proportional_richness <- sweep(MDiv, 2, max_richness, FUN = "/", check.margin = FALSE)
MDiv$Multidiversity <- rowMeans(proportional_richness, na.rm = TRUE)

###****** Calculate EMF index
###*https://cran.r-project.org/web/packages/multifunc/vignettes/multifunc_biodepth.html
EMF <- EMF %>% mutate(logCH4 = -1*logCH4 +max(logCH4, na.rm=T)) #revert the totC variable
EMF <- EMF %>% mutate(logCHLa = -1*logCHLa +max(logCHLa, na.rm=T)) #revert the totC variable
EMF_vars <- c("logCH4","logCHLa","logBMzp","logBMmi","Plants")
EMF <- cbind(EMF, getStdAndMeanFunctions(EMF, EMF_vars)) #scale your data by dividing by the max of each column 
EMF <- EMF[,c(6:11)]
EMF_MDiv = cbind(MDiv,EMF)
names(EMF_MDiv)[c(5,11)]=paste(c("MDiv","EMF"))

#Transform your factors
EF_SEM$logArea <- log(EF_SEM$Area)
EF_SEM$logDepth <- log(EF_SEM$Depth)
EF_SEM$logTN <- log(EF_SEM$TN)
EF_SEM$logTP <- log(EF_SEM$TP)
#scale the variables you cannot transform. 
EF_SEM_vars<-c("Cropland", "Longitude", "Latitude","MAT","AI", "Aridity")
EF_SEM<-cbind(EF_SEM, getStdAndMeanFunctions(EF_SEM, EF_SEM_vars)) #scale your data by dividing by the max of each column 

#reorder your variables, combine the database and add AI as binomial variable
TF <- cbind(EMF_MDiv,EF_SEM)
TF <- TF[,c(1:11,35:44, 13, 32)]
TF$Ard <- ifelse(EF_SEM$AI >= 0.650000, 0, "1")

# Save your database
write.csv(TF, file = "FullDB_T.csv")


#####################################################################################******
################################# Full set of functions #############################******
#####################################################################################******
#set up your tables
EF_SEM <- read_excel("EMF_rawData.xlsx", sheet = "FullDB_AllF")

EF <- EF_SEM #N=105
EF = EF %>% drop_na(S,k,SedTraps)
PondID_Country <- EF[ , c(1,2)]
EF <- column_to_rownames(EF, var="Pond_ID")
EF <- EF[ , -c(1,23,24)]

###****** Permutate the NA of your funtions and taxonomic groups
EMF = EF[,c(1:8)]
EMF <- mice(EMF, method = "pmm", m = 5)  # Permutation with m = 5. Increase if necessary
EMF <- complete(EMF)  # extract the data
#Transform your EF
EMF$logS <- log(EMF$S)
EMF$logk <- sqrt(EMF$k)
EMF$logST <- log(EMF$SedTraps)
EMF$logCH4 <- sqrt(EMF$CH4)
EMF$logCHLa <- log(EMF$CHLa)
EMF$logBMzp <- log(EMF$BMzp)
EMF$logBMmi <- log(EMF$BMmi)
EMF = EMF[,c(9:15,8)]

###****** Calculate multidiversity
###*This is calculated in the same way that weight EMF from Wang 2019
MDiv = EF[,c(9:12)]
max_richness <- apply(MDiv, 2, max, na.rm = TRUE)
proportional_richness <- sweep(MDiv, 2, max_richness, FUN = "/", check.margin = FALSE)
MDiv$Multidiversity <- rowMeans(proportional_richness, na.rm = TRUE)

###****** Calculate EMF index
###*https://cran.r-project.org/web/packages/multifunc/vignettes/multifunc_biodepth.html
EMF <- EMF %>% mutate(logCH4 = -1*logCH4 +max(logCH4, na.rm=T)) #revert the totC variable
EMF <- EMF %>% mutate(logCHLa = -1*logCHLa +max(logCHLa, na.rm=T)) #revert the totC variable
EMF_vars<-c("logS","logk","logST","logCH4","logCHLa","logBMzp","logBMmi","Plants")
EMF<-cbind(EMF, getStdAndMeanFunctions(EMF, EMF_vars)) #scale your data by dividing by the max of each column 
EMF_MDiv = cbind(MDiv,EMF)
EMF_MDiv = EMF_MDiv[,c(1:5,14:22)]
names(EMF_MDiv)[c(5,14)]=paste(c("MDiv","EMF"))

#Transform your factors
EF_SEM$logArea <- log(EF_SEM$Area)
EF_SEM$logDepth <- log(EF_SEM$Depth)
EF_SEM$logTN <- log(EF_SEM$TN)
EF_SEM$logTP <- log(EF_SEM$TP)
#scale the variables you cannot transform. 
EF_SEM_vars<-c("Cropland", "Longitude", "Latitude","MAT", "AI", "Aridity")
EF_SEM<-cbind(EF_SEM, getStdAndMeanFunctions(EF_SEM, EF_SEM_vars)) #scale your data by dividing by the max of each column 
EF_SEM = EF_SEM %>% drop_na(S,k,SedTraps)
#reorder your variables, combine the database and add AI as binomial variable
TF <- cbind(EMF_MDiv,EF_SEM)
TF <- TF[,c(1:14,41:50,16,38)]
TF$Ard <- ifelse(EF_SEM$AI >= 0.650000, 0, "1")
# Save your database
write.csv(TF, file = "FullDB_AllF_T.csv")