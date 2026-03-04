setwd("~/Dropbox/PONDERFUL/Database/data/Ecosystem functions/Analisis/All/SEM/Article/R1_11Sep2025/R1_documents/To Submit/DataScript")
source("~/Dropbox/PONDERFUL/Database/data/Ecosystem functions/Analisis/All/SEM/Article/Documents to submit/Database_Script/PackagesFuntions.R")

############################################################################## Main figures
#******Fig 1 - Aridity Index Map

#******Fig 2a, c-f - B-EMF relationship
#Upload you table 
TF <- read_excel("EMF_rawDataT.xlsx",sheet = "TF")
#create the model and check the results
model<-lme(EMF ~ BioMI, #substitute MDiv by the different diversities (MDiv, BioMP, BioZP, BioMI, BioAmph)
             random = list(~ Country|Pondscape, ~ 1|Year), #random = list(~ Country|Pondscape, ~ 1|Year), #(~ 1|Country, ~ Month|Year),
             control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
             weights = varIdent(form = ~ 1|Year), #Specifies that the residual variance differs across different months within each year. This allows for modeling heteroscedasticity.
             na.action = na.omit,
             data= TF)
summary <- papeR::prettify(summary(model), signif.stars=F, digits=3, smallest.pval = 0.001)
summary
r.squaredGLMM(model, null.RE = TRUE)
#Check model assumptions
plot(model)
qqnorm(resid(model))
#Generate overall fixed predicted values for lme models
TF$multidiversity_EF_fixed = predict(model, newdata = TF, level = 0)
##Extracting Confidence intervals for lme model (fixed effect) for plotting
des_multidiversity_EF = model.matrix(formula(model)[-2], TF)
predvar_multidiversity_EF = diag(des_multidiversity_EF %*% vcov(model) %*% t(des_multidiversity_EF))
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
Corr <- cor(TF[,c(7,3:6)],TF[,c(13,8:12)], method="pearson") 
corrplot(Corr,method="color", addCoef.col = "black",col.lim = c(-0.5, 0.5),is.corr = FALSE,cl.pos = "b")

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
ANCmodel <- aov(EMF ~ MDiv*as.character(Ard), data=TF)
summary(ANCmodel)
#model assumptions
qqnorm(resid(ANCmodel), main="Normal Q-Q Plot")

#plot the the model with arity as a categorical factor
TF$Ard=as.character(TF$Ard)
#Remplace MDiv by BioMP, BioZP, BioMI and BioAmph
model<-lme(EMF ~ BioAmph*Ard,
           random = list(~ Country|Pondscape), 
           control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
           na.action = na.omit,
           data= TF)
summary (model)
r.squaredGLMM(model, null.RE = TRUE)
#plotting. Change the biodiversity in the model 
TF$multidiversity_MEF_fixed <- predict(model, newdata = TF, level = 0)
des_multidiversity_MEF <- model.matrix(~ BioAmph * as.character(Ard), data = TF)
predvar_multidiversity_MEF <- diag(des_multidiversity_MEF %*% vcov(model) %*% t(des_multidiversity_MEF))
TF$lower_multidiversity_MEF <- TF$multidiversity_MEF_fixed - 2 * sqrt(predvar_multidiversity_MEF)
TF$upper_multidiversity_MEF <- TF$multidiversity_MEF_fixed + 2 * sqrt(predvar_multidiversity_MEF)
#Remplace MDiv by BioMP, BioZP, BioMI and BioAmph
ggplot(TF, aes(x = BioAmph, y = EMF)) +
  geom_point(aes(color = Ard), size = 1.5, alpha = 0.6) +
  geom_line(aes(y = multidiversity_MEF_fixed, color = Ard, group = Ard), linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower_multidiversity_MEF, ymax = upper_multidiversity_MEF, fill = Ard, group = Ard), alpha = 0.1) +
  theme_minimal() +
  labs(y = "EMF", color = "Aridity", fill = "Aridity") +
  scale_color_manual(values = c("0" = "#99008C", "1" = "#F6A245")) +
  scale_fill_manual(values = c("0" = "#99008C", "1" = "#F6A245"))

#dividing the database
TFard <- filter(TF, Ard=="1")
TFhum <- filter(TF, Ard=="0")
  
model<-lme(EMF ~ BioMP,
             random = list(~ Country|Pondscape), 
             control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
             na.action = na.omit,
             data= TFard)
summary (model)
r.squaredGLMM(model, null.RE = TRUE)

#******Fig 4. SEM analysis
myModel <- ' 
   EMF ~ MDiv + FishPA + logTN + logTP + logArea + logDepth + MAT.std + Aridity.std + Cropland.std + Longitude.std + Latitude.std
   MDiv ~ FishPA + logTN + logTP + logArea + logDepth + MAT.std + Aridity.std + Cropland.std + Longitude.std + Latitude.std
   FishPA ~ logTN + logTP + logArea + logDepth + MAT.std + Aridity.std + Longitude.std + Latitude.std
   logTN ~ Cropland.std + logArea + logDepth + MAT.std + Aridity.std + Longitude.std + Latitude.std
   logTP ~ Cropland.std + logArea + logDepth + MAT.std + Aridity.std + Longitude.std + Latitude.std
   logArea ~ Cropland.std + MAT.std + Aridity.std + Longitude.std + Latitude.std 
   logDepth ~ Cropland.std + MAT.std + Aridity.std + Longitude.std + Latitude.std
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


###################################################################### Supplementary figures
#******Sup Fig.1. Supporting figure for the methane and Chl-a negative relationship with ZP diversity
#upload row data
EF_SEM <- read_excel("EMF_rawData.xlsx", sheet = "FullDB")

#transform variables for normality assumptions
EF_SEM$logCHLa <- log(EF_SEM$CHLa)
EF_SEM$sqrCH4 <- sqrt(EF_SEM$CH4) #EF_SEM$logTotC <- log(EF_SEM$TotC)
EF_SEM$logBMzp <- log(EF_SEM$BMzp)
EF_SEM$logBMmi <- log(EF_SEM$BMmi)
EF_SEM$logPlants <- log(EF_SEM$Plants)
EF_SEM$logTP <- log(EF_SEM$TP)
EF_SEM$logTN <- log(EF_SEM$TN)

#create the correlations between variables
SupFig1 = EF_SEM[,c(9,26,25,24,30,29)]
ggpairs(SupFig1, 
        lower = list(continuous = wrap(lowerFn, method = "lm")),
        diag = list("continuous" = function(data, mapping, ...){ggally_text(rlang::as_label(mapping$x),col="black") + theme_void()}),
        upper = list(continuous = wrap("cor", method="pearson",color = "black")))
#create the models
#Model1 with zooplankton richness as dependent variable of biomass, methane and chlorophyll-a
Date <- read_excel("EMF_rawDataT.xlsx", sheet = "TF")
Date <- Date[ , c(26:28)]
EF_SEM <- cbind(EF_SEM,Date)
model1<-lme(BioZP ~ logBMzp+sqrCH4+logCHLa, #substitute by the different diversities (MP, ZP, MI, Amph)
            random = list(~ 1|Country, ~ Month|Year), 
            control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
            weights = varIdent(form = ~ Month|Year), #Specifies that the residual variance differs across different months within each year. This allows for modeling heteroscedasticity.
            na.action = na.omit,
            data= EF_SEM)
summary1 <- papeR::prettify(summary(model1), signif.stars=F, digits=3, smallest.pval = 0.001)
summary1
r.squaredGLMM(model1, null.RE = TRUE)
#Model2 with chlorophyll-a as dependn variable of methane, nitrogen and phosphorous concentration. 
model2<-lme(logCHLa ~ sqrCH4+logTP+logTN, #substitute by the different diversities (MP, ZP, MI, Amph)
            random = list(~ 1|Country, ~ Month|Year), #random = list(~ 1|Country, ~ Month|Year),
            control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
            weights = varIdent(form = ~ Month|Year), #Specifies that the residual variance differs across different months within each year. This allows for modeling heteroscedasticity.
            na.action = na.omit,
            data= EF_SEM)
summary2 <- papeR::prettify(summary(model2), signif.stars=F, digits=3, smallest.pval = 0.001)
summary2
r.squaredGLMM(model1, null.RE = TRUE)

#******Sup Fig.2. Conceptual figure of threshold analysis

#******Supplementary Fig 3. Threshold analysis of biodiversity 
PercentCol <- c("90%" = "#BAE9C6", "75%" = "#7DD2BF", "50%" = "#37ADB6", "25%" = "#0078A2", "10%" ="#2D3184")
EMF_MDiv <- read_csv("EMF_MDiv.csv")
biodiversity_list <- c(MDiv = "MDiv",MP = "BioMP",ZP = "BioZP",MI = "BioMI",Amph = "BioAmph")# Biodiversity variable names and labels for plotting
allVars <- c("logCH4.std", "logCHLa.std", "logBMzp.std", "logBMmi.std", "Plants.std")# Environmental function variables to include
target_thresholds <- c(0.1, 0.25, 0.5, 0.75, 0.9)# Target thresholds to extract (with tolerance handling)
tolerance <- 1e-6
ThrEFfinal <- data.frame()

# Loop through each biodiversity metric
for (taxon in names(biodiversity_list)) {
  biodiv_col <- biodiversity_list[[taxon]]
  DB <- EMF_MDiv[, c("Pond_ID", biodiv_col, allVars)]# Subset relevant columns: Pond_ID, current biodiversity, and EMF vars
  names(DB)[1:2] <- c("Pond_ID", "Diversity")# Rename for compatibility
  vars <- whichVars(DB, allVars)# Get variable indices for EMF vars
  # Calculate thresholds
  DBThresh <- getFuncsMaxed(DB, vars, threshmin = 0.05, threshmax = 0.99, prepend = c("Pond_ID", "Diversity"), maxN = 20)
  DBThresh$percent <- paste0(100 * DBThresh$thresholds, "%")# Add percent label
  # Filter selected thresholds using numeric tolerance
  gcPlot <- DBThresh[sapply(DBThresh$thresholds, function(x) any(abs(x - target_thresholds) < tolerance)), ]
  gcPlot$Taxa <- taxon# Add taxon name
  # Combine into final data frame
  ThrEFfinal <- bind_rows(ThrEFfinal, gcPlot)
}
write.csv(ThrEFfinal, file = "ThrEFfinal.csv",row.names = FALSE)
#Create your plot
ggplot(ThrEFfinal, aes(x = Diversity, y = funcMaxed, color = factor(percent))) +
  stat_smooth(method = "lm", se = FALSE, method.args = list(family = quasipoisson(link = "identity"))) + 
  #geom_jitter(width = 0.1, height = 0.1, alpha = 0.2, size=0.5) + 
  scale_color_manual(name="percent", values=PercentCol) +
  ylab(expression("Number of Functions >= Threshold")) +
  theme(legend.position="none") + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  xlab("Diversity") + theme_classic() + 
  facet_wrap(~Taxa, scales = "free")

#Statistical analysis
DBThresh<-filter(ThrEFfinal, Taxa=="MI")#Substitute by Mdiv, MP, ZP, MI and Amph
modelDB<-glm(funcMaxed ~ Diversity, data=subset(DBThresh, DBThresh$thresholds=="0.1"), family=quasipoisson(link="identity"))
summary(modelDB)
modelDB<-glm(funcMaxed ~ Diversity, data=subset(DBThresh, DBThresh$thresholds=="0.25"), family=quasipoisson(link="identity"))
summary(modelDB)
modelDB<-glm(funcMaxed ~ Diversity, data=subset(DBThresh, DBThresh$thresholds=="0.5"), family=quasipoisson(link="identity"))
summary(modelDB)
modelDB<-glm(funcMaxed ~ Diversity, data=subset(DBThresh, DBThresh$thresholds=="0.75"), family=quasipoisson(link="identity"))
summary(modelDB)
modelDB<-glm(funcMaxed ~ Diversity, data=subset(DBThresh, DBThresh$thresholds=="0.9"), family=quasipoisson(link="identity"))
summary(modelDB)
#modelDB<-lm(funcMaxed ~ Diversity, data=subset(DBThresh, DBThresh$thresholds=="0.9"))
#summary(modelDB)



#******Sup Fig.4. Boxplot with diversities and functions. 
#create the long format table with the specific variables
SupFig4=TF[,c(1,3:13,25)]
SupFig4 = SupFig4 %>% pivot_longer(!c(Pond_ID,Ard), names_to = "Variable", values_to = "Values")
#anova analysis
Model <- aov(MDiv ~ as.character(Ard), data=TF)
summary(Model)
#plot the resutls
ggplot(SupFig4, aes(x=as.character(Ard), y=Values)) + #label = spp
  geom_boxplot(position=position_dodge(1), outlier.shape = NA, color="black", fill="grey", alpha=0.2) +
  geom_jitter(width = 0.3, alpha = 0.5, aes(colour = as.character(Ard))) +
  theme(legend.position="none") + facet_wrap(~Variable, scales = "free")

#******Supplementary Figure 5. Standardized total effect from SEM 
#Values extracted from the results file
STF <- data.frame(
  Variable = c("Latitude", "Longitude","TN","TP","Croplants","Fish P/A","Depth","Area","Aridity","MAT","Multidiversity"),
  StandTotEffect = c(-0.11,0.005,-0.037,-0.031,-0.034,-0.048,-0.024,0.009,0.050,0.17,0.278)
)
STF$Variable <- factor(STF$Variable, levels = STF$Variable)
ggplot(STF, aes(x = Variable, y = StandTotEffect)) + 
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() + labs(y = "Standardized total effects") + theme_minimal()

#******Sup Fig.6. Pearson correlations functions and biodiversity
#CorrAll <- TF[,c(14,9:13,8,4:7)] #div-8,4:7 or Fun-14,9:13
CorrAll <- TF[,c(14,9:13,25,24,23,15,16,21,20,22,19,18)] #div-8,4:7 or Fun-14,9:13
CorrMatrix <- corr_coef(CorrAll, method = "pearson") #change by pearson or spearman
plot(CorrMatrix, type = "upper", reorder = FALSE, diag = TRUE) +
  scale_fill_gradient2(low = "#B2182B", high ="#2166AC", midpoint = 0,
                       limits = c(-0.8,0.8), space = "Lab",
                       name="Pearson\nCorrelation")

#******Sup Fig.7. Relationship between richness and shannon index 
SI <- read_excel("ShannonIndex.xlsx")
ggscatter(SI, x = "ShanZP", y = "BioZP", size = 2,
          add.params = list(position = position_jitter(w = .1, h = .1),color = "orange", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", font.label = c(20, "plain"))
ggscatter(SI, x = "ShanMI", y = "BioMI", size = 2,
          add.params = list(position = position_jitter(w = .1, h = .1),color = "orange", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", font.label = c(20, "plain"))

#******Sup Fig.8. Pearson correlation
CorrMDiv <- TF[,c(3:6)] 
CorrMatrix <- corr_coef(CorrMDiv, method = "pearson") #change by pearson or spearman
plot(CorrMatrix, type = "upper", reorder = FALSE, diag = TRUE) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high ="#2166AC", midpoint = 0,
                       limits = c(-0.5,0.5), space = "Lab",
                       name="Pearson\nCorrelation")
CorrEF <- EF_SEM[,c(24:28)]
CorrEF <- CorrEF[complete.cases(CorrEF), ]
CorrMatrix <- corr_coef(CorrEF, method = "pearson") #change by pearson or spearman
plot(CorrMatrix, type = "upper", reorder = FALSE, diag = TRUE) +
  scale_fill_gradient2(low = "#B2182B", high ="#2166AC", midpoint = 0,
                       limits = c(-0.5,0.5), space = "Lab",
                       name="Pearson\nCorrelation")

#******Supplementary Figure 9. SEM with all the function and a data subset 
TF <- read_csv("TF.csv")
TF <- TF[,c(1,12)]
TF_all <- read_csv("TF_allF.csv")
TF_all <- TF_all[,c(1,15)]
TFcorr <- merge(TF_all, TF, by = c("...1"))
ggscatter(TFcorr, x = "EMF.x", y = "EMF.y", size = 2,
          add.params = list(position = position_jitter(w = .1, h = .1),color = "orange", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", font.label = c(20, "plain"))

#******Supplementary Figure 10. SEM with all the funciton and a data subset
TF <- read_excel("EMF_rawDataT.xlsx",sheet = "TF_subsetDB")
myModel <- ' 
   EMF ~ MDiv + FishPA + logTN + logTP + logArea + logDepth + MAT.std + Aridity.std + Cropland.std + Longitude.std + Latitude.std
   MDiv ~ FishPA + logTN + logTP + logArea + logDepth + MAT.std + Aridity.std + Cropland.std + Longitude.std + Latitude.std
   FishPA ~ logTN + logTP + logArea + logDepth + MAT.std + Aridity.std + Longitude.std + Latitude.std
   logTN ~ Cropland.std + logArea + logDepth + MAT.std + Aridity.std + Longitude.std + Latitude.std
   logTP ~ Cropland.std + logArea + logDepth + MAT.std + Aridity.std + Longitude.std + Latitude.std
   logArea ~ Cropland.std + MAT.std + Aridity.std + Longitude.std + Latitude.std 
   logDepth ~ Cropland.std + MAT.std + Aridity.std + Longitude.std + Latitude.std
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
parameterEstimates(fit, se = TRUE, zstat = TRUE, pvalue = TRUE, ci = TRUE, standardized = FALSE, rsquare = TRUE)
# get the test statistic for the model fit
fit_FD <- cfa(myModel, data = TF, se="none")

#******Supplementary Table 7. Bootstrapping for all the possible combinations of functions
# Create a vector of variables
#IMPORTANT!!! Run from the Pre-cleaning analysis until the line #33
EMF_vars <- c("logCH4", "logCHLa", "logBMzp", "logBMmi", "Plants")
all_results <- list()
i <- 1

for (k in 1:length(EMF_vars)) {
  combs <- combn(EMF_vars, k, simplify = FALSE)
  for (vars in combs) {
    result <- getStdAndMeanFunctions(EMF, vars)
    # Convert to data frame and safely add the combination label
    result_df <- as.data.frame(result)
    result_df$Combination <- paste(vars, collapse = ", ")
    # Store in list
    all_results[[i]] <- result_df
    i <- i + 1
  }
}

final_df <- bind_rows(all_results)
#extract the EMF index calculated with all the combinations of the individual functions:
final_df$PondNumb<-c(rep(1:206,31))
final_df <- final_df[,c(8,3,2)]
BT_EMF = final_df %>% pivot_wider(names_from = Combination, values_from = meanFunction)
#remove those combinations with only one functions. The final database have 26 combinations
BT_EMF <- BT_EMF[,c(7:32)]
#merge the database with the MDiv
BT_EMF = cbind(MDiv,BT_EMF)
Date <- read_excel("EMF_rawDataT.xlsx", sheet = "TF")
Date <- Date[ , c(26:27,2)]
BT_EMF = cbind(Date,BT_EMF) #add year and month to include as random effect in the models
BT_EMF <- BT_EMF[,c(1:3,8,9:34)]
#to be sure of the calculations the column with the combination of the five functions need to be the same results than the original analysis
names(BT_EMF)[c(4,30)]=paste(c("MDiv","EMF"))
#remove spaces from the names of the columns
names(BT_EMF) <- gsub(", ", "_", names(BT_EMF))

#create the model to check the the relationship between MDiv and the different combination of functions
# List of response variable names (columns in BT_EMF)
response_vars <- c("logCH4_logCHLa","logCH4_logBMzp","logCH4_logBMmi","logCH4_Plants",
  "logCHLa_logBMzp","logCHLa_logBMmi","logCHLa_Plants","logBMzp_logBMmi",
  "logBMzp_Plants","logBMmi_Plants","logCH4_logCHLa_logBMzp", "logCH4_logCHLa_logBMmi",
  "logCH4_logCHLa_Plants","logCH4_logBMzp_logBMmi","logCH4_logBMzp_Plants","logCH4_logBMmi_Plants",
  "logCHLa_logBMzp_logBMmi","logCHLa_logBMzp_Plants","logCHLa_logBMmi_Plants","logBMzp_logBMmi_Plants",
  "logCH4_logCHLa_logBMzp_logBMmi","logCH4_logCHLa_logBMzp_Plants",
  "logCH4_logCHLa_logBMmi_Plants","logCH4_logBMzp_logBMmi_Plants", "logCHLa_logBMzp_logBMmi_Plants","EMF")
# Initialize list to store model results
model_results <- list()

for (resp in response_vars) {
  cat("Fitting model for response:", resp, "\n")
  # Create formula: <response> ~ MDiv
  model_formula <- as.formula(paste(resp, "~ MDiv"))
  # Try fitting the model
  model_try <- try({
    model <- lme(
      fixed = model_formula,
      random = list(~1 | Country, ~Month | Year),
      control = list(msMaxIter = 200, opt = "optim", msVerbose = TRUE),
      weights = varIdent(form = ~ Month | Year),
      na.action = na.omit,
      data = BT_EMF)
    #Summarise the model and extraxt the R2
    summary_pretty <- papeR::prettify(summary(model), signif.stars = FALSE, digits = 3, smallest.pval = 0.001)
    r2 <- r.squaredGLMM(model, null.RE = TRUE)
    #create the list of results you want to extract
    list(response = resp,model = model,summary = summary_pretty,r_squared = r2)
  }, silent = TRUE)
  # Store result
  model_results[[resp]] <- if (inherits(model_try, "try-error")) {
    list(error = model_try)
  } else {
    model_try
  }
}

#to create a table with the results you want to extract from the model:
# Create an empty list to store rows
model_summary_list <- list()

for (resp in names(model_results)) {
  res <- model_results[[resp]]
  if (!is.null(res$error)) next #Skip failed models
  coefs <- summary(res$model)$tTable #Extract coefficients and p-values
  r2 <- res$r_squared #Extract R²
  # Add a row to the summary list
  model_summary_list[[resp]] <- data.frame(
    Response = resp,
    Estimate = coefs["MDiv", "Value"],
    Std_Error = coefs["MDiv", "Std.Error"],
    DF = coefs["MDiv", "DF"],
    t_value = coefs["MDiv", "t-value"],
    p_value = coefs["MDiv", "p-value"],
    Marginal_R2 = r2[1],
    Conditional_R2 = r2[2],
    row.names = NULL
  )
}

# Combine into a single dataframe and save the results 
model_summary_df <- do.call(rbind, model_summary_list)
write.csv(model_summary_df, file = "BT_EMFresults.csv", row.names = FALSE)



