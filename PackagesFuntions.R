# packages ----------------------------------------------------------------
library(readr)
library(readxl)
library(dplyr)
library(data.table)
library(ggpubr)
library(metan)
library(GGally)
library(ggplot2)
library(corrplot)
library(tidyverse)
library(reshape2)
library(viridis) 
library(grid)
library(gridExtra)
library(lme4)
library(car)
library(nlme)
library(MuMIn)
library(multifunc)
library(mice)
library(ggpubr)
library(jtools)
library(lavaan)

lowerFn <- function(data, mapping, method = "lm", ...) { p <- ggplot(data = data, mapping = mapping) +
  geom_jitter(width = 0.3, alpha = 0.5) +
  geom_smooth(method = method, color = "black")
p}
