setwd("~/Dropbox/PONDERFUL/Database/data/Ecosystem functions/Analisis/All/SEM/Article/R1_11Sep2025/NewAnalysis/R_scripts")
source("~/Dropbox/PONDERFUL/Database/data/Ecosystem functions/Analisis/All/SEM/Article/Documents to submit/Database_Script/PackagesFuntions.R")

#set up your tables
EF_SEM <- read_excel("EMF_rawData_DKSP.xlsx")

EF <- EF_SEM #N=60
PondID_Country <- EF[ , c(1,2,3)] #extract pond id, country and pondscapes
EF <- column_to_rownames(EF, var="Pond_ID")
EF <- EF[ , -c(1,2,21,22)] #extract the binomial variables

###****** Permutate the NA of your functions and taxonomic group 
EMF = EF[,c(1:5)]
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

###****** Calculate multidiversity
#This is calculated in the same way that weight EMF from Wang 2019
MDiv = EF[,c(6:9)]
max_richness <- apply(MDiv, 2, max, na.rm = TRUE)
proportional_richness <- sweep(MDiv, 2, max_richness, FUN = "/", check.margin = FALSE)
MDiv$Multidiversity <- rowMeans(proportional_richness, na.rm = TRUE)

###****** Calculate EMF index
#Corr <- cor(EMF_MDiv, method="pearson") 
#corrplot(Corr,method="color", addCoef.col = "black",col.lim = c(-1, 1),is.corr = FALSE,cl.pos = "b")
###*https://cran.r-project.org/web/packages/multifunc/vignettes/multifunc_biodepth.html
EMF <- EMF %>% mutate(logCH4 = -1*logCH4 +max(logCH4, na.rm=T)) #revert the totC variable
EMF <- EMF %>% mutate(logCHLa = -1*logCHLa +max(logCHLa, na.rm=T)) #revert the CHLa variable
#EMF <- EMF %>% mutate(NEP = -1*NEP +max(NEP, na.rm=T)) #revert the NEP variable
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
TF <- TF[,c(12:14,1:11,43:45,36:40,31:34,25)]
TF$Ard <- ifelse(TF$AI >= 0.650000, 0, "1")

# Save your database
write.csv(TF, file = "FullDB_T_DKSP.csv", )



setwd("~/Dropbox/PONDERFUL/Database/data/Ecosystem functions/Analisis/All/SEM/Article/R1_11Sep2025/NewAnalysis/R_scripts")
source("~/Dropbox/PONDERFUL/Database/data/Ecosystem functions/Analisis/All/SEM/Article/Documents to submit/Database_Script/PackagesFuntions.R")

############################################################################## Main figures
#******Fig 1 - Aridity Index Map

#******Fig 2a, c-f - B-EMF relationships
#Upload you table 
#TF <- read_excel("EMF_rawDataT.xlsx",sheet = "TF")
TF <- read_csv("FullDB_T_DKSP.csv")
TF=TF[,c(2:28)]
#create the model and check the results
model<-lme(EMF ~ MDiv, #substitute MDiv by the different diversities (bioMP, bioZP, bioMI, bioAmph)
           random = list(~ 1|Pondscape, ~ 1|Year), #random = list(~ Pondscape|Country, ~ Month|Year),
           control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
           weights = varIdent(form = ~ 1|Year), #Specifies that the residual variance differs across different months within each year. This allows for modeling heteroscedasticity.
           na.action = na.omit,
           data= TF)
summary <- papeR::prettify(summary(model), signif.stars=F, digits=3, smallest.pval = 0.001)
summary
r.squaredGLMM(model, null.RE = TRUE)
#Check model assumptions
#plot(model) xx
#qqnorm(resid(model)) xx
#Generate overall fixed predicted values for lme models
TF$multidiversity_EF_fixed = predict(model, newdata = TF, level = 0)
##Extracting Confidence intervals for lme model (fixed effect) for plotting
des_multidiversity_EF = model.matrix(formula(model)[-2], TF)
predvar_multidiversity_EF = diag(des_multidiversity_EF %*% vcov(model) %*%
                                   t(des_multidiversity_EF) )
TF$lower_multidiversity_EF = with(TF,multidiversity_EF_fixed - 2*sqrt(predvar_multidiversity_EF))
TF$upper_multidiversity_EF = with(TF,multidiversity_EF_fixed + 2*sqrt(predvar_multidiversity_EF))
#created the figure from the extracted predicted values 
#substitute MDiv by the different diversities (bioMP, bioZP, bioMI, bioAmph)
ggplot(TF, aes(x = MDiv, y = EMF)) +
  geom_point(size = 1.5, position = position_jitter(width = 0, height = 0.00), stroke = 0.3, alpha = 0.6) +
  geom_line(aes(y = multidiversity_EF_fixed), color = "orange", size = 1) +
  geom_ribbon(aes(y = NULL, ymin = lower_multidiversity_EF, ymax = upper_multidiversity_EF), alpha = 0.08) +
  theme_minimal()

#******Fig 2b - Pearson correlations functions and biodiversity
Corr <- cor(TF[,c(8,4:7)],TF[,c(14,9:13)], method="pearson") 
corrplot(Corr,method="color", addCoef.col = "black",col.lim = c(-0.5, 0.5),is.corr = FALSE,cl.pos = "b")
ggscatter(EF_SEM, x = "BioZP", y = "logCH4", size = 2,
          add.params = list(position = position_jitter(w = .1, h = .1),color = "orange", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", font.label = c(20, "plain"))

#******Fig 3a - B-EMF relationship across a gradient of aridity - Standardized effect
#scale all your variables
data_model<-TF %>% mutate(
  MDiv= (MDiv - min(MDiv, na.rm = TRUE))/max(MDiv - min(MDiv, na.rm = TRUE)),
  BioMP= (BioMP - min(BioMP, na.rm = TRUE))/max(BioMP - min(BioMP, na.rm = TRUE)),
  BioZP= (BioZP - min(BioZP, na.rm = TRUE))/max(BioZP - min(BioZP, na.rm = TRUE)),
  BioMI= (BioMI - min(BioMI, na.rm = TRUE))/max(BioMI - min(BioMI, na.rm = TRUE)),
  BioAmph= (BioAmph - min(BioAmph, na.rm = TRUE))/max(BioAmph - min(BioAmph, na.rm = TRUE)))
#create the models
#model by using multidiversity
model<-lme(EMF ~ MDiv*AI.std,
           random = list(~ 1|Pondscape), 
           control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
           na.action = na.omit,
           data= data_model)
summary <- papeR::prettify(summary(model), signif.stars=F, digits=3, smallest.pval = 0.001)
colnames(summary)[2] <- "Estimate"
summary
r.squaredGLMM(model, null.RE = TRUE)
#model1 by using macrophytes richness
model1<-lme(EMF ~ BioMP*AI.std,
            random = list(~ 1|Pondscape), 
            control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
            na.action = na.omit,
            data= data_model)
summary1 <- papeR::prettify(summary(model1), signif.stars=F, digits=3, smallest.pval = 0.001)
colnames(summary1)[2] <- "Estimate"
summary1
r.squaredGLMM(model1, null.RE = TRUE)
#model2 by using zooplankton richness
model2<-lme(EMF ~ BioZP*AI.std,
            random = list(~ 1|Pondscape), 
            control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
            na.action = na.omit,
            data= data_model)
summary2 <- papeR::prettify(summary(model2), signif.stars=F, digits=3, smallest.pval = 0.001)
colnames(summary2)[2] <- "Estimate"
summary2
r.squaredGLMM(model2, null.RE = TRUE)
#model3 by using macroinvertebrates richness
model3<-lme(EMF ~ BioMI*AI.std,
            random = list(~ 1|Pondscape), 
            control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
            na.action = na.omit,
            data= data_model)
summary3 <- papeR::prettify(summary(model3), signif.stars=F, digits=3, smallest.pval = 0.001)
colnames(summary2)[2] <- "Estimate"
summary3
r.squaredGLMM(model3, null.RE = TRUE)
#model4 by using amphibians richness
model4<-lme(EMF ~ BioAmph*AI.std,
            random = list(~ 1|Pondscape), 
            control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
            na.action = na.omit,
            data= data_model)
summary4 <- papeR::prettify(summary(model4), signif.stars=F, digits=3, smallest.pval = 0.001)
colnames(summary2)[2] <- "Estimate"
summary4
r.squaredGLMM(model4, null.RE = TRUE)
#create the figure with the standardized fixed effects
plot_summs(model,model1,model2,model3,model4,
           scale = TRUE,
           robust = TRUE,
           ci_level = 0.9,
           colors = "Qual1",
           inner_ci_level = NULL,
           plot.distributions = F,
           rescale.distributions = T,
           omit.coefs = c("(Intercept)", "Arity.std"),
           model.names = NULL)

#******Fig 3b-f, Comparison between Arid and humid sites
#Ancova. Remplace MDiv by BioMP, BioZP, BioMI and BioAmph
ANCmodel <- aov(EMF ~ BioAmph*as.character(Ard), data=TF)
summary(ANCmodel)
#model assumptions
#qqnorm(resid(ANCmodel), main="Normal Q-Q Plot")
#plot the the model with arity as a categorical factor
TF$Ard=as.character(TF$Ard)
#Remplace MDiv by BioMP, BioZP, BioMI and BioAmph
model<-lme(EMF ~ MDiv*Ard,
           random = list(~ 1|Pondscape), 
           control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
           na.action = na.omit,
           data= TF)
summary(model)
TF$multidiversity_MEF_fixed <- predict(model, newdata = TF, level = 0)
des_multidiversity_MEF <- model.matrix(~ MDiv * as.character(Ard), data = TF)
predvar_multidiversity_MEF <- diag(des_multidiversity_MEF %*% vcov(model) %*% t(des_multidiversity_MEF))
TF$lower_multidiversity_MEF <- TF$multidiversity_MEF_fixed - 2 * sqrt(predvar_multidiversity_MEF)
TF$upper_multidiversity_MEF <- TF$multidiversity_MEF_fixed + 2 * sqrt(predvar_multidiversity_MEF)
#Remplace MDiv by BioMP, BioZP, BioMI and BioAmph
ggplot(TF, aes(x = MDiv, y = EMF)) +
  geom_point(aes(color = Ard), size = 1.5, alpha = 0.6) +
  geom_line(aes(y = multidiversity_MEF_fixed, color = Ard, group = Ard), linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower_multidiversity_MEF, ymax = upper_multidiversity_MEF, fill = Ard, group = Ard), alpha = 0.1) +
  theme_minimal() + theme(legend.position="none") + 
  labs(x = "Multidiversity", y = "EMF", color = "Ard", fill = "Ard") +
  scale_color_manual(values = c("0" = "#99008C", "1" = "#F6A245")) +
  scale_fill_manual(values = c("0" = "#99008C", "1" = "#F6A245"))


#modelos individuales
DK=filter(TF, Country=="Denmark")
SP=filter(TF, Country=="Spain")
model<-lme(EMF ~ MDiv,
            random = list(~ 1|Pondscape), 
            control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
            na.action = na.omit,
            data= SP)
summary <- papeR::prettify(summary(model), signif.stars=F, digits=3, smallest.pval = 0.001)
colnames(summary2)[2] <- "Estimate"
summary
r.squaredGLMM(model, null.RE = TRUE)

ggscatter(DK, x = "BioZP", y = "EMF", size = 2,
          add.params = list(position = position_jitter(w = .1, h = .1),color = "orange", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", font.label = c(20, "plain"))
