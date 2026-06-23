#analysis of the data
TF <- read.csv("~/Dropbox/PONDERFUL/Database/data/Ecosystem functions/Analisis/All/SEM/EMF_T_DKSP.csv")
names(TF)[c(1)]=paste(c("Pond_ID"))
#Add month and year for the models
EMF_rawDataT <- read_excel("DB_Script/EMF_rawDataT.xlsx")
EMF_rawDataT <- EMF_rawDataT %>%
  filter(Country %in% c("Denmark", "Spain")) %>%
  select(1,tail(names(EMF_rawDataT), 3))
TF <- merge(TF, EMF_rawDataT, by = c("Pond_ID"))

#******Fig 2a, c-f - B-EMF relationship across a gradient of aridity
#create the model and check the results
model<-lme(EMF ~ BioMP, #substitute MDiv by the different diversities (bioMP, bioZP, bioMI, bioAmph)
           random = list(~ 1|Year), #random = list(~ 1|Country, ~ Month|Year),
           control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
           weights = varIdent(form = ~ 1|Year), #Specifies that the residual variance differs across different months within each year. This allows for modeling heteroscedasticity.
           na.action = na.omit,
           data= TF)
model <- lm(EMF ~ BioAmph, data=TF)
summary <- papeR::prettify(summary(model), signif.stars=F, digits=3, smallest.pval = 0.001)
summary
r.squaredGLMM(model, null.RE = TRUE)
#Check model assumptions
plot(model)
qqnorm(resid(model))
#Generate overall fixed predicted values for lme models
TF$multidiversity_EF_fixed <- predict(model, newdata = TF, level = 0)
# Correct variance for CI
des_multidiversity_EF <- model.matrix(formula(model)[-2], TF)
predvar_multidiversity_EF <- diag(des_multidiversity_EF %*% vcov(model) %*% 
                                    t(des_multidiversity_EF))
TF$lower_multidiversity_EF <- TF$multidiversity_EF_fixed - 2*sqrt(predvar_multidiversity_EF)
TF$upper_multidiversity_EF <- TF$multidiversity_EF_fixed + 2*sqrt(predvar_multidiversity_EF)
# Sort before plotting
TF <- TF[order(TF$MDiv), ]
#created the figure from the extracted predicted values 
#substitute MDiv by the different diversities (bioMP, bioZP, bioMI, bioAmph)
ggplot(TF, aes(x = BioAmph, y = EMF)) +
  geom_point(size = 1.5, position = position_jitter(width = 0, height = 0.00), stroke = 0.3, alpha = 0.6) +
  geom_line(aes(y = multidiversity_EF_fixed), color = "orange", size = 1) +
  geom_ribbon(aes(y = NULL, ymin = lower_multidiversity_EF, ymax = upper_multidiversity_EF), alpha = 0.1) +
  theme_minimal()

#******Fig 2b - Pearson correlations functions
Corr <- cor(TF[,c(6,2:5)],TF[,c(12,7:11)], method="pearson") 
corrplot(Corr,method="color", addCoef.col = "black",col.lim = c(-0.6, 0.6),is.corr = FALSE,cl.pos = "b")

library(metan)
library(GGally)
colcheck = TF[,c(12,7:11,5,4,3,2,6)] #Diversity [], functions [14:19] colcheck = M2[,c(1,38,20:21,23:28)]
#colcheck = M5[,c("tM5","FishPA","ECELS.std","MeanT.std","Aridity.std","logArea","logDepth","logTN","logTP","Croplant.std")]
CorrMatrix <- corr_coef(colcheck, method = "pearson") #change by pearson or spearman
plot(CorrMatrix,type = "upper",reorder = FALSE, diag = TRUE, legend.title = "Pearson's\nCorrelation")

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
           random = list(~ 1|SamplingDate), 
           control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
           na.action = na.omit,
           data= data_model)
summary <- papeR::prettify(summary(model), signif.stars=F, digits=3, smallest.pval = 0.001)
colnames(summary)[2] <- "Estimate"
summary
r.squaredGLMM(model, null.RE = TRUE)
#model1 by using macrophytes richness
model1<-lme(EMF ~ BioMP*AI.std,
            random = list(~ 1|SamplingDate), 
            control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
            na.action = na.omit,
            data= data_model)
summary1 <- papeR::prettify(summary(model1), signif.stars=F, digits=3, smallest.pval = 0.001)
colnames(summary1)[2] <- "Estimate"
summary1
r.squaredGLMM(model1, null.RE = TRUE)
#model2 by using zooplankton richness
model2<-lme(EMF ~ BioZP*AI.std,
            random = list(~ 1|SamplingDate), 
            control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
            na.action = na.omit,
            data= data_model)
summary2 <- papeR::prettify(summary(model2), signif.stars=F, digits=3, smallest.pval = 0.001)
colnames(summary2)[2] <- "Estimate"
summary2
r.squaredGLMM(model2, null.RE = TRUE)
#model3 by using macroinvertebrates richness
model3<-lme(EMF ~ BioMI*AI.std,
            random = list(~ 1|SamplingDate), 
            control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
            na.action = na.omit,
            data= data_model)
summary3 <- papeR::prettify(summary(model3), signif.stars=F, digits=3, smallest.pval = 0.001)
colnames(summary2)[2] <- "Estimate"
summary3
r.squaredGLMM(model3, null.RE = TRUE)
#model2 by using amphibians richness
model4<-lme(EMF ~ BioAmph*AI.std,
            random = list(~ 1|SamplingDate), 
            control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
            na.action = na.omit,
            data= data_model)
summary4 <- papeR::prettify(summary(model4), signif.stars=F, digits=3, smallest.pval = 0.001)
colnames(summary2)[2] <- "Estimate"
summary4
r.squaredGLMM(model4, null.RE = TRUE)
#create the figure with the standarised fixed effects
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
ANCmodel <- aov(EMF ~ BioMP*as.character(Ard), data=TF)
summary(ANCmodel)
#model assumptions
qqnorm(resid(ANCmodel), main="Normal Q-Q Plot")
#plot the the model with arity as a categorical factor
data_model$Ard=as.character(data_model$Ard)
TF$Ard=as.character(TF$Ard)
#Remplace MDiv by BioMP, BioZP, BioMI and BioAmph
model<-lme(EMF ~ MDiv*Ard,
           random = list(~ 1|SamplingDate), 
           control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
           na.action = na.omit,
           data= TF)
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
  theme_minimal() +
  labs(y = "EMF", color = "Aridity", fill = "Aridity") +
  scale_color_manual(values = c("0" = "#99008C", "1" = "#F6A245")) +
  scale_fill_manual(values = c("0" = "#99008C", "1" = "#F6A245"))

TFdry <- TF[TF$Ard %in% c("1"), ]
TFhum <- TF[TF$Ard %in% c("0"), ]
model<-lme(EMF ~ BioMP,
           random = list(~ 1|SamplingDate), 
           control=list(msMaxIter=200,opt = "optim",msVerbose=TRUE),
           na.action = na.omit,
           data= TFdry)
summary <- papeR::prettify(summary(model), signif.stars=F, digits=3, smallest.pval = 0.001)
summary
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

