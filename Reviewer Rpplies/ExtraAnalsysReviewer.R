source("~/Dropbox/PONDERFUL/Database/data/Ecosystem functions/Analisis/All/SEM/Article/Documents to submit/Database_Script/PackagesFuntions.R")

####################### 1st reviewer
####################### PCA for land use
#######################
###****** Extract PCA axes 
EFtOUTL <- read_excel("LUvisual.xlsx", sheet = "PCA_check")
LU <- column_to_rownames(EFtOUTL, var="Pond_ID")
## scree plot
LuPCA <- PCA(LU) #Change the database (LU5, LU100, LU500)
#create a visual plot
fviz_screeplot(LuPCA, ncp = 10) + theme_classic()
LuPCA <- prcomp(LU, center = TRUE, scale. = TRUE) #Change by LU5, LU100 and LU500
## variable plot
fviz_pca_var(LuPCA, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07", "#FF0000"),
             repel = TRUE,   # avoid text overlapping
             axes = c(1, 2)) # choose PCs to plot

####################### 1st reviewer
####################### SP-DK comparison, it is in independent scripts. 
#######################

####################### 2nd reviewer
####################### new Pearson correlation between the functions and the biodiversity, without revert the functions.
####################### Fig. 4 MAIN TEXT
#This script is part of the DB_Preparation database with all data together
#set up your tables
setwd("~/Dropbox/PONDERFUL/Database/data/Ecosystem functions/Analisis/All/SEM/Article/Documents to submit/Database_Script")
EF_SEM <- read_excel("EMF_rawData.xlsx", sheet = "FullDB")
EF <- EF_SEM #N=206
PondID_Country <- EF[ , c(1,2)] #extract pond id and country
EF <- column_to_rownames(EF, var="Pond_ID")
EF <- EF[ , -c(1,20,21)] #extract the binomial variables
###****** Permutate the NA of your functions
EMF = EF[,c(1:5)]
EMF <- mice(EMF, method = "pmm", m = 5)  # Permutation with m = 5. Increase if necessary
EMF <- complete(EMF)  # extract the data
#Transform your EF to reach normality
EMF$logCH4 <- log(EMF$CH4) 
EMF$logCHLa <- log(EMF$CHLa)
EMF$logBMzp <- log(EMF$BMzp)
EMF$logBMmi <- log(EMF$BMmi)
EMF$KPlogCH4 <- log(EMF$CH4) 
EMF$KPlogCHLa <- log(EMF$CHLa)
EMF$KPlogBMzp <- log(EMF$BMzp)
EMF$KPlogBMmi <- log(EMF$BMmi)
EMF = EMF[,c(6:9,5,10:13)]
###****** Calculate multidiversity
#This is calculated in the same way that weight EMF from Wang 2019
MDiv = EF[,c(6:9)]
max_richness <- apply(MDiv, 2, max, na.rm = TRUE)
proportional_richness <- sweep(MDiv, 2, max_richness, FUN = "/", check.margin = FALSE)
MDiv$Multidiversity <- rowMeans(proportional_richness, na.rm = TRUE)
###****** Calculate EMF index
###*https://cran.r-project.org/web/packages/multifunc/vignettes/multifunc_biodepth.html
EMF <- EMF %>% mutate(logCH4 = -1*logCH4 +max(logCH4, na.rm=T)) #revert the CH4 variable
EMF <- EMF %>% mutate(logCHLa = -1*logCHLa +max(logCHLa, na.rm=T)) #revert the Chl-a variable
EMF_vars <- c("logCH4","logCHLa","logBMzp","logBMmi","Plants")
EMF <- cbind(EMF, getStdAndMeanFunctions(EMF, EMF_vars)) #scale your data by dividing by the max of each column 
EMF_MDiv = cbind(MDiv,EMF)
names(EMF_MDiv)[c(5,20)]=paste(c("MDiv","EMF"))

#******Fig 2b - Pearson correlations functions and biodiversity
Corr <- cor(EMF_MDiv[,c(5,1:4)],EMF_MDiv[,c(20,11:19)], method="pearson") 
corrplot(Corr,method="color", addCoef.col = "black",col.lim = c(-0.5, 0.5),is.corr = FALSE,cl.pos = "b")

colcheck = EMF_MDiv[,c(20,11:19,4,3,2,1,5)] 
CorrMatrix <- corr_coef(colcheck, method = "pearson") #change by pearson or spearman
plot(CorrMatrix,type = "upper",reorder = FALSE, diag = TRUE, legend.title = "Pearson's\nCorrelation")

####################### 3rd reviewer
####################### Pondscape and country B-EMF relationship 
#######################
setwd("~/Dropbox/PONDERFUL/Database/data/Ecosystem functions/Analisis/All/SEM/Article/Documents to submit/Database_Script")
TF <- read_excel("EMF_rawDataT.xlsx",sheet = "TF")

#Ancova. Remplace MDiv by BioMP, BioZP, BioMI and BioAmph
ANCmodel <- aov(EMF ~ MDiv*Country, data=TF)
ANCmodel <- lm(EMF ~ MDiv+Country, data=TF)
summary(ANCmodel)
#model assumptions
qqnorm(resid(ANCmodel), main="Normal Q-Q Plot")

#General model with all countries included. [Remplace MDiv by BioMP, BioZP, BioMI and BioAmph]
model<-lm(EMF ~ MDiv*Country,
          na.action = na.omit,
          data= TF)
TF$multidiversity_MEF_fixed <- predict(model, newdata = TF, level = 0)
des_multidiversity_MEF <- model.matrix(~ MDiv * Country, data = TF)
predvar_multidiversity_MEF <- diag(des_multidiversity_MEF %*% vcov(model) %*% t(des_multidiversity_MEF))
TF$lower_multidiversity_MEF <- TF$multidiversity_MEF_fixed - 2 * sqrt(predvar_multidiversity_MEF)
TF$upper_multidiversity_MEF <- TF$multidiversity_MEF_fixed + 2 * sqrt(predvar_multidiversity_MEF)
#Remplace MDiv by BioMP, BioZP, BioMI and BioAmph
ggplot(TF, aes(x = MDiv, y = EMF)) +
  geom_point(aes(color = Country), size = 1.5, alpha = 0.6) +
  geom_line(aes(y = multidiversity_MEF_fixed, color = Country, group = Country), linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower_multidiversity_MEF, ymax = upper_multidiversity_MEF, fill = Country, group = Country), alpha = 0.1) +
  theme_minimal() + xlim(0.1, 0.7) + ylim(0.1, 0.7) + theme(legend.position="none")

#change of spatial scale for individual countries and pondscape as variable [Remplace MDiv by BioMP, BioZP, BioMI and BioAmph]
XXCountry<-filter(TF, Country=="Switzerland") #change by the country. This is a selection only per country
#Each color represent a country
PondscapeCol <- c("#AE123A", "#CC243C", "#EB4F48","#FA8D76", "#FDAB9B") #Belgium
PondscapeCol <- c("#E5803F", "#E9A328", "#E6D04F","#E2E6BD", "#E3E400") #Denmark
PondscapeCol <- c("#09622A", "#176F31", "#519E53","#81B274") #Germany
PondscapeCol <- c("#038692", "#009380", "#32AE7C","#5DBE75", "#94D268") #Spain
PondscapeCol <- c("#0019FF", "#005DFF", "#0080FF","#00A2FF", "#00B9BC","#00E5FF") #Switzerland
PondscapeCol <- c("#3E0689","#663A96", "#88499D", "#A65CA4", "#C3A2DE", "#F1D5FA") #Turkey
PondscapeCol <- c("#A42175", "#B72678", "#EA4374","#E2438A", "#ED567D") #UK

model<-lm(EMF ~ MDiv*Pondscape, data= XXCountry)
summary(model)
XXCountry$multidiversity_MEF_fixed <- predict(model, newdata = XXCountry, level = 0)
des_multidiversity_MEF <- model.matrix(~ MDiv*Pondscape, data = XXCountry)
predvar_multidiversity_MEF <- diag(des_multidiversity_MEF %*% vcov(model) %*% t(des_multidiversity_MEF))
XXCountry$lower_multidiversity_MEF <- XXCountry$multidiversity_MEF_fixed - 2 * sqrt(predvar_multidiversity_MEF)
XXCountry$upper_multidiversity_MEF <- XXCountry$multidiversity_MEF_fixed + 2 * sqrt(predvar_multidiversity_MEF)
ggplot(XXCountry, aes(x = MDiv, y = EMF)) +
  geom_point(aes(color = Pondscape), size = 1.5, alpha = 0.6) +
  geom_line(aes(y = multidiversity_MEF_fixed, color = Pondscape, group = Pondscape), linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower_multidiversity_MEF, ymax = upper_multidiversity_MEF, group = Pondscape), alpha = 0.05) +
  theme_minimal() + scale_color_manual(name="Pondscape", values=PondscapeCol) + xlim(0.1, 0.7) + ylim(0.1, 0.7) + theme(legend.position="none")

#to perform the individual models per landscape and check the correlation. Most of the correlations are not signification
#because the number of data are too small
TF$Pondscape <- factor(TF$Pondscape, levels = c("Antwerp","Flemish_Brabant","Limburg_2","Limburg_1","East_Flanders","Laven","Avernako",
                                                "Aero","Holstebro","RyHule","Muncheberg","Schoneiche","Lietzen","Quillow","Albera",
                                                "Osona","Selva","Garrtoxa","Gavarres","Jussy","Rhone","Versoix","Champagne","Meyrin",
                                                "Seymaz","Alpagut","Imrendi","Akbas","Ayas Yolu","Sorgun", "Karacaoren",
                                                "Norfolk Bodham-Baconsthorpe","Norfolk Horningtoft-Brisley","Norfolk Heydon","Lancashire",
                                                "Cheshire"))
XXPondsp<-filter(TF, Pondscape=="Flemish_Brabant")
model<-lm(EMF ~ MDiv*Pondscape, data= TF)
summary(model)

####################### 3rd reviewer
####################### Weighted VS Unweighted multifunctionality
#######################
TF <- read_excel("~/Dropbox/PONDERFUL/Database/data/Ecosystem functions/Analisis/All/SEM/Article/Documents to submit/Database_Script/EMF_rawDataT.xlsx")
functions<-TF  %>% dplyr::select("logCH4.std","logCHLa.std","logBMzp.std","logBMmi.std","Plants.std")
function_nrs <- c(1:5) # Those functions analyzed for ecosystem function MF
names(functions)[function_nrs]
#DENDOGRAMA
functions_complete = complete(functions)
functions[,function_nrs] <- functions_complete
functions_matrix <- t(as.matrix(functions[,function_nrs]))
functions_matrix <- scale(functions_matrix)
d <- dist(functions_matrix, method = "euclidean")
dendrogram <- hclust(d, method = "complete" )
plot(dendrogram, cex = 0.6, hang = -1) 

function_nrs2 <- function_nrs
names(functions)[function_nrs2]
loading_values <- c(0.5,0.33,0.33,0.33,0.5)

#Threshold level of 'high' functioning
threshold <- 0.5
#rescale your function
mf_data_scaled <- data.matrix(functions)
for (i in seq_along(function_nrs2)) {
  x <- functions[[ function_nrs2[i] ]]
  maximum <- mean(sort(x, decreasing = TRUE)[1:5])
  #mf_data_scaled[,i] <- rep(0,length(mf_data_scaled[,i]))
  high_values <- which(x >= threshold * maximum)
  mf_data_scaled[high_values, i] <- loading_values[i]
}
#Calculate weight ecosystem multifunctionality: sum of loadings above threshold, divided by 3 (maximum possible score), to get a value between 0 and 1
functions$weighted_multifunctionality <- rowSums(mf_data_scaled) / 2

MF = TF[,c(1,13)]
df <- bind_cols(functions, MF)
library("ggpubr")
ggscatter(df, x = "weighted_multifunctionality", y = "EMF", size = 2,
          add.params = list(position = position_jitter(w = .1, h = .1),color = "orange", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", font.label = c(20, "plain")) +
          xlab("Weighted EMF") + ylab("Unweighted EMF")

####################### 2nd reviewer
####################### Hydroperiod addition
#######################
setwd("~/Dropbox/PONDERFUL/Database/data/Ecosystem functions/Analisis/All/SEM/Article/Documents to submit/Database_Script")
TF <- read_excel("EMF_rawDataT.xlsx",sheet = "TF")
EF_SEM <- read_excel("EMF_rawData.xlsx", sheet = "FullDB")
MF = EF_SEM[,c(1,22)]
TF <- bind_cols(TF, MF)

CorrAll <- TF[,c(14,15,17:24,31)]
CorrAll <- EF_SEM[,c(12:22)]
CorrMatrix <- corr_coef(CorrAll, method = "pearson") #change by pearson or spearman
plot(CorrMatrix, type = "upper", reorder = FALSE, diag = TRUE) +
  scale_fill_gradient2(low = "#B2182B", high ="#2166AC", midpoint = 0,
                       limits = c(-0.8,0.8), space = "Lab",
                       name="Pearson\nCorrelation")

myModel <- ' 
   EMF ~ MDiv + FishPA + logTN + logTP + logArea + logDepth + MAT.std + Aridity.std + Cropland.std + Longitude.std + Latitude.std + Hydrop
   MDiv ~ FishPA + logTN + logTP + logArea + logDepth + MAT.std + Aridity.std + Cropland.std + Longitude.std + Latitude.std + Hydrop
   FishPA ~ logTN + logTP + logArea + logDepth + MAT.std + Aridity.std + Longitude.std + Latitude.std + Hydrop
   logTN ~ Cropland.std + logArea + logDepth + MAT.std + Aridity.std + Longitude.std + Latitude.std + Hydrop
   logTP ~ Cropland.std + logArea + logDepth + MAT.std + Aridity.std + Longitude.std + Latitude.std + Hydrop
   logArea ~ Cropland.std + MAT.std + Aridity.std + Longitude.std + Latitude.std + Hydrop
   logDepth ~ Cropland.std + MAT.std + Aridity.std + Longitude.std + Latitude.std + Hydrop
   Hydrop ~ Cropland.std + logArea + logDepth + MAT.std + Aridity.std
   Cropland.std ~ MAT.std + Aridity.std + Longitude.std + Latitude.std
   Aridity.std ~ Longitude.std + Latitude.std
   MAT.std ~ Longitude.std + Latitude.std 
   
#Residual correlations Between environmental variables
   logTN ~~ logTP
   logArea ~~ logDepth
   Longitude.std ~~ Latitude.std 
   MAT.std ~~ Aridity.std
'
###Analyzing model output
fit <- sem(model = myModel, data = TF, fixed.x = F) 
summary(fit, standardize = TRUE, rsq = TRUE, fit.measures = TRUE)
# get the test statistic for the model fit
fit_FD <- cfa(myModel, data = TF, se="none")
summary(fit_FD, fit.measures = TRUE)
# get the test statistic for the original sample
T.orig_FD <- fitMeasures(fit, "CFI")
# bootstrap to get bootstrap test statistics
T.boot_FD <- bootstrapLavaan(fit_FD, R=1000, type="bollen.stine",
                             FUN=fitMeasures, fit.measures="CFI") #Between 5 to 10% of failure it is acceptable

T.boot_FD_clean <- na.omit(T.boot_FD)
# compute a bootstrap based p-value
pvalue.boot_FD <- length(which(T.boot_FD_clean < T.orig_FD))/length(T.boot_FD_clean)
1-pvalue.boot_FD 
# extract the parameters
parameterEstimates(fit, se = TRUE, zstat = TRUE, pvalue = TRUE, ci = TRUE, standardized = FALSE, rsquare = TRUE)

####################### reviewer 1
####################### trade off. figure sup mat 8
#######################
CorrMDiv <- TF[,c(3:6)] 
CorrMatrix <- corr_coef(CorrMDiv, method = "pearson") #change by pearson or spearman
plot(CorrMatrix, type = "upper", reorder = FALSE, diag = TRUE) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high ="#2166AC", midpoint = 0,
                       limits = c(-0.5,0.5), space = "Lab",
                       name="Pearson\nCorrelation")

EF_SEM <- read_excel("EMF_rawData.xlsx", sheet = "FullDB")
#transform variables for normality assumptions
EF_SEM$logCHLa <- log(EF_SEM$CHLa)
EF_SEM$sqrCH4 <- sqrt(EF_SEM$CH4) #EF_SEM$logTotC <- log(EF_SEM$TotC)
EF_SEM$logBMzp <- log(EF_SEM$BMzp)
EF_SEM$logBMmi <- log(EF_SEM$BMmi)
EF_SEM$logPlants <- log(EF_SEM$Plants)
EF_SEM$logTP <- log(EF_SEM$TP)
EF_SEM$logTN <- log(EF_SEM$TN)
CorrEF <- EF_SEM[,c(24:28)]
CorrEF <- CorrEF[complete.cases(CorrEF), ]
CorrMatrix <- corr_coef(CorrEF, method = "pearson") #change by pearson or spearman
plot(CorrMatrix, type = "upper", reorder = FALSE, diag = TRUE) +
  scale_fill_gradient2(low = "#B2182B", high ="#2166AC", midpoint = 0,
                       limits = c(-0.5,0.5), space = "Lab",
                       name="Pearson\nCorrelation")

#upload row data
EF_SEMall <- read_excel("EMF_rawData.xlsx", sheet = "FullDB_AllF")
#transform variables for normality assumptions
EF_SEMall$logCHLa <- log(EF_SEMall$CHLa)
EF_SEMall$sqrCH4 <- sqrt(EF_SEMall$CH4) #EF_SEM$logTotC <- log(EF_SEM$TotC)
EF_SEMall$logBMzp <- log(EF_SEMall$BMzp)
EF_SEMall$logBMmi <- log(EF_SEMall$BMmi)
EF_SEMall$logPlants <- log(EF_SEMall$Plants)
EF_SEMall$logS <- log(EF_SEMall$S)
EF_SEMall$logk <- log(EF_SEMall$k)
EF_SEMall$logSedTraps <- log(EF_SEMall$SedTraps)

CorrEF <- EF_SEMall[,c(28:35)]
CorrEF <- CorrEF[complete.cases(CorrEF), ]
CorrMatrix <- corr_coef(CorrEF, method = "pearson") #change by pearson or spearman
plot(CorrMatrix, type = "upper", reorder = FALSE, diag = TRUE) +
  scale_fill_gradient2(low = "#B2182B", high ="#2166AC", midpoint = 0,
                       limits = c(-0.5,0.5), space = "Lab",
                       name="Pearson\nCorrelation")

####################### 3rd reviewer
####################### Remove year 2022 and provide and interactions
#######################
TF <- read_excel("EMF_rawDataT.xlsx",sheet = "TF")
model<-lme(EMF ~ MDiv * Year, #substitute MDiv by the different diversities (BioMP, BioZP, BioMI, BioAmph)
           random = list(~ Country|Pondscape), #Also tested: random = list(~ Country|Pondscape, ~ Month|Year),
           control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
          # weights = varIdent(form = ~ 1|Year), #Specifies that the residual variance differs across different months within each year. This allows for modeling heteroscedasticity.
           na.action = na.omit,
           data= TF)
summary <- papeR::prettify(summary(model), signif.stars=F, digits=3, smallest.pval = 0.001)
summary
r.squaredGLMM(model, null.RE = TRUE)
#Check model assumptions
plot(model)
qqnorm(resid(model))
TF$multidiversity_EF_fixed = predict(model, newdata = TF, level = 0)
des_multidiversity_EF = model.matrix(~ MDiv * Year, data = TF)
predvar_multidiversity_EF = diag(des_multidiversity_EF %*% vcov(model) %*% t(des_multidiversity_EF) )
TF$lower_multidiversity_EF = with(TF,multidiversity_EF_fixed - 2*sqrt(predvar_multidiversity_EF))
TF$upper_multidiversity_EF = with(TF,multidiversity_EF_fixed + 2*sqrt(predvar_multidiversity_EF))
ggplot(TF, aes(x = MDiv, y = EMF)) +
  geom_point(aes(color = Year), size = 1.5, alpha = 0.6) +
  geom_line(aes(y = multidiversity_EF_fixed, color = Year, group = Year), linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower_multidiversity_EF, ymax = upper_multidiversity_EF, fill = Year, group = Year), alpha = 0.1) +
  theme_minimal() + theme(legend.position="none")

TF21 = filter(TF, Year=="2021")
model<-lme(EMF ~ MDiv, #substitute MDiv by the different diversities (BioMP, BioZP, BioMI, BioAmph)
           random = list(~ Country|Pondscape), #Also tested: random = list(~ Country|Pondscape, ~ Month|Year),
           control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
           # weights = varIdent(form = ~ 1|Year), #Specifies that the residual variance differs across different months within each year. This allows for modeling heteroscedasticity.
           na.action = na.omit,
           data= TF21)
summary <- papeR::prettify(summary(model), signif.stars=F, digits=3, smallest.pval = 0.001)
summary
r.squaredGLMM(model, null.RE = TRUE)
#Check model assumptions
plot(model)
qqnorm(resid(model))

TF21$multidiversity_EF_fixed = predict(model, newdata = TF21, level = 0)
des_multidiversity_EF = model.matrix(formula(model)[-2], TF21)
predvar_multidiversity_EF = diag(des_multidiversity_EF %*% vcov(model) %*% t(des_multidiversity_EF) )
TF21$lower_multidiversity_EF = with(TF21,multidiversity_EF_fixed - 2*sqrt(predvar_multidiversity_EF))
TF21$upper_multidiversity_EF = with(TF21,multidiversity_EF_fixed + 2*sqrt(predvar_multidiversity_EF))
ggplot(TF21, aes(x = MDiv, y = EMF)) +
  geom_point(size = 1.5, position = position_jitter(width = 0, height = 0.00), stroke = 0.3, alpha = 0.6) +
  geom_line(aes(y = multidiversity_EF_fixed), color = "blue", size = 1) +
  geom_ribbon(aes(y = NULL, ymin = lower_multidiversity_EF, ymax = upper_multidiversity_EF), alpha = 0.08) +
  theme_minimal()

####################### 3rd reviewer
####################### To normalize variables - Z-Score
#######################
#set up your tables
EF_SEM <- read_excel("EMF_rawData.xlsx")
EF <- EF_SEM
###****** Permutate the NA of your functions and taxonomic group 
EMF = EF[,c(4:8)]
EMF <- mice(EMF, method = "pmm", m = 5)  # Permutation with m = 5. Increase if necessary
EMF <- complete(EMF)  # extract the data
#Transform your EF to reach normality
#ggdensity(EMF$NEP) #ggqqplot(EMF$NEP) #shapiro.test(EMF$Straps)
EMF$logCH4 <- log(EMF$CH4) 
EMF$logCHLa <- log(EMF$CHLa)
EMF$logBMzp <- log(EMF$BMzp)
EMF$logBMmi <- log(EMF$BMmi)
#EMF$logStraps <- log(EMF$Straps) #EMF = EMF[,c(8:11,5,12,7)]
EMF = EMF[,c(6:9,5)]
###****** Calculate EMF index
EMF <- EMF %>% mutate(logCH4 = -1*logCH4 +max(logCH4, na.rm=T)) #revert the totC variable
EMF <- EMF %>% mutate(logCHLa = -1*logCHLa +max(logCHLa, na.rm=T)) #revert the CHLa variable
EMF_vars <- c("logCH4","logCHLa","logBMzp","logBMmi","Plants")
EMF <- cbind(EMF, getStdAndMeanFunctions(EMF, EMF_vars)) #scale your data by dividing by the max of each column 
EMFminmax <- EMF[,c(6:11)]


###****** Permutate the NA of your functions and taxonomic group 
EMF = EF[,c(4:8)]
EMF <- mice(EMF, method = "pmm", m = 5)  # Permutation with m = 5. Increase if necessary
variables_interes <- complete(EMF)  # extract the data
#created the transformation quantil-quantil function
transformar_columna <- function(x) {
  perc <- rank(x, na.last = "keep") / (sum(!is.na(x)) + 1)  #keep NA
  result <- ifelse(!is.na(x), qnorm(perc), NA)  
  return(result)
}
#applied the function "transformar_columna" to your selected variables and extract the table
###****** Calculate EMF index
variables_interesN <- apply(variables_interes, 2, transformar_columna)
EMFzscore <- as.data.frame(variables_interesN)
EMFzscore <- EMFzscore %>% mutate(CH4 = -1*CH4 +max(CH4, na.rm=T)) #revert the totC variable
EMFzscore <- EMFzscore %>% mutate(CHLa = -1*CHLa +max(CHLa, na.rm=T)) #revert the CHLa variable
EMFzscore$EMFzscore <- rowMeans(EMFzscore[, 1:5]) #scale your data by dividing by the max of each column 

##Check the correlation between. 
EMF2meth = cbind(EMFzscore,EMFminmax)
write.csv(EMF2meth, file = "EMF2meth.csv", )
library(metan)
library(GGally)
CorrMatrix <- corr_coef(EMF2meth, method = "pearson") #change by pearson or spearman
plot(CorrMatrix,type = "upper",reorder = FALSE, diag = TRUE, legend.title = "Pearson's\nCorrelation")
ggscatter(EMF2meth, x = "EMFzscore", y = "meanFunction", size = 2,
          add.params = list(position = position_jitter(w = .1, h = .1),color = "orange", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", font.label = c(20, "plain")) +
          xlab("Z-score EMF") + ylab("MinMax EMF")

#check the different models with the two scaling methods. 
EMF2methRev3 <- read_excel("EMF2methRev3.xlsx")
#create the model and check the results
model<-lme(EMFzscore ~ MDiv, #substitute MDiv by the different diversities (bioMP, bioZP, bioMI, bioAmph)
           random = list(~ 1|Pondscape, ~ 1|Year), #random = list(~ Pondscape|Country, ~ Month|Year),
           control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
           weights = varIdent(form = ~ 1|Year), #Specifies that the residual variance differs across different months within each year. This allows for modeling heteroscedasticity.
           na.action = na.omit,
           data= EMF2methRev3)
summary <- papeR::prettify(summary(model), signif.stars=F, digits=3, smallest.pval = 0.001)
summary
r.squaredGLMM(model, null.RE = TRUE)
#Check model assumptions
#plot(model) xx
#qqnorm(resid(model)) xx
#Generate overall fixed predicted values for lme models
EMF2methRev3$multidiversity_EF_fixed = predict(model, newdata = EMF2methRev3, level = 0)
##Extracting Confidence intervals for lme model (fixed effect) for plotting
des_multidiversity_EF = model.matrix(formula(model)[-2], EMF2methRev3)
predvar_multidiversity_EF = diag(des_multidiversity_EF %*% vcov(model) %*%
                                   t(des_multidiversity_EF) )
EMF2methRev3$lower_multidiversity_EF = with(EMF2methRev3,multidiversity_EF_fixed - 2*sqrt(predvar_multidiversity_EF))
EMF2methRev3$upper_multidiversity_EF = with(EMF2methRev3,multidiversity_EF_fixed + 2*sqrt(predvar_multidiversity_EF))
#created the figure from the extracted predicted values 
#substitute MDiv by the different diversities (bioMP, bioZP, bioMI, bioAmph)
ggplot(EMF2methRev3, aes(x = MDiv, y = EMFzscore)) +
  geom_point(size = 1.5, position = position_jitter(width = 0, height = 0.00), stroke = 0.3, alpha = 0.6) +
  geom_line(aes(y = multidiversity_EF_fixed), color = "orange", size = 1) +
  geom_ribbon(aes(y = NULL, ymin = lower_multidiversity_EF, ymax = upper_multidiversity_EF), alpha = 0.08) +
  theme_minimal()


####################### 3rd reviewer
####################### Threshold analysis
#######################
#upload the data
EMF_MDiv <- read_csv("EMF_MDiv.csv")
allVars <- c("logCH4.std", "logCHLa.std", "logBMzp.std", "logBMmi.std", "Plants.std")
#the FuncMax IS THE MAXIMAL FUNCION THAT ARE IN EACH THRESHOLD AND IS CREATED BY THE NEXT FUNCTION: 
DBThresh<-getFuncsMaxed(EMF_MDiv, allVars, threshmin=0.05, threshmax=0.99, prepend=c("Pond_ID"), maxN=20)
#Given that variation, let’s look at the entire spread of thresholds.
DBThresh$percent <- 100*DBThresh$thresholds
#change the biodivesity by MDiv, BioMP, BioZP, BioMI and BioAmph
ggplot(data=DBThresh, aes(x=BioAmph, y=funcMaxed, group=percent)) +
  ylab(expression("Number of Functions" >= Threshold)) +
  xlab("Species Richness") +
  stat_smooth(method="glm", 
              method.args = list(family=quasipoisson(link="log")), 
              lwd=0.8, fill=NA, aes(color=percent)) +
  theme_bw(base_size=14) +
  scale_color_gradient(name="Percent of \nMaximum", low="blue", high="red")

#Standard desviation
LinearSlopes<-getCoefTab(funcMaxed ~ BioAmph,data = DBThresh, coefVar = "BioAmph")
# Plot the values of the diversity slope at different levels of the threshold
Slopes <- ggplot(LinearSlopes, aes(x=thresholds*100, y = estimate, ymax = estimate + 1.96 * std.error, ymin = estimate - 1.96 * std.error)) +
  geom_ribbon(fill="grey50") + geom_point() +
  ylab("Change in Number of Functions per Addition of 1 Species\n") + xlab("\nThreshold (%)") +
  geom_abline(intercept=0, slope=0, lwd=1, linetype=2) + theme_bw(base_size=14)
Slopes

LinearSlopes<-getCoefTab(funcMaxed ~ BioAmph,
                                data = DBThresh, 
                                coefVar = "BioAmph")
IDX <- getIndices(LinearSlopes, DBThresh, funcMaxed ~ BioAmph)
IDX

