# ---
#title: "Formatting UK Biobank data for the refreshed data"
#author: "Jean Shin"
#date: "2024-04-09"
# 
# loading libraries and define functions --------------------------------------

rm(list=ls())
setwd("~/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP")

library(tidyverse)
library(data.table)
library(patchwork)
library(ggcorrplot)
library(ggpubr)
library(EnvStats)

inormal <- function(x)
{
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}
remove.outliers_grubbs <- function(dat, varnames){
  require(outliers)
  count.na0 <- count.na1 <- c()
  
  for(varname in varnames){
    x = dat[[varname]]
    count.na0 = c(count.na0,sum(is.na(x)))
    keep.going <- T
    while(keep.going){
      test <- grubbs.test(x,opposite = T)
      print(test)
      if(test$p.value<0.05){
        if(str_detect(test$alternative,"lowest")){
          cat("-------------------------------------------------------------------")
          cat(paste('\nReplacing the lowest value of \'',varname,'\' with NA.\n',sep=''))
          cat("-------------------------------------------------------------------\n")
          x[which.min(x)] <- NA
        }else if(str_detect(test$alternative,"highest")){
          cat("-------------------------------------------------------------------")
          cat(paste('\nReplacing the highest value of \'',varname,'\' with NA.\n',sep=''))
          cat("-------------------------------------------------------------------")
          x[which.max(x)] <- NA
        }
      }
      test2 <- grubbs.test(x,opposite = F)
      print(test2)
      if(test2$p.value<0.05){
        if(str_detect(test2$alternative,"lowest")){
          cat("-------------------------------------------------------------------")
          cat(paste('\nReplacing the lowest value of \'',varname,'\' with NA.\n',sep=''))
          cat("-------------------------------------------------------------------")
          x[which.min(x)] <- NA
        }else if(str_detect(test2$alternative,"highest")){
          cat("-------------------------------------------------------------------")
          cat(paste('\nReplacing the highest value of \'',varname,'\' with NA.\n',sep=''))
          cat("-------------------------------------------------------------------")
          x[which.max(x)] <- NA
        }
      }
      
      keep.going=test$p.value<0.05 | test2$p.value<0.05
    }
    count.na1 = c(count.na1,sum(is.na(x)))
    dat[[varname]] <- x
  }
  counts.NA = data.frame(var=varnames,before=count.na0,after=count.na1)
  ret = list(counts.NA = counts.NA, cleaed_data=dat)
  ret
}
  
# d_brain + ICV ---------------------------------------------------------------
d_brain_ICV = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb/data/ukb677610_d_brain_unrelated.txt')

# type 2 diabetes -------------------------------------------------------------
d_T2D = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb/data/SummaryDiabetesDiagnosis_baseline.txt')

# bmi and covariates ----------------------------------------------------------
d_bmi_covs = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb/data/ukb677610_d_bmi_covariates.txt')

# apoE and geno PCs -----------------------------------------------------------
d_geno = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb/data/ukb677610_d_apoE4_genoPC.txt')

# brain pheno names -----------------------------------------------------------
library(RCurl)
BrainPhenoNames <- getURL("https://raw.githubusercontent.com/jshinb/neuroCHARGE_scripts/main/Obesity_CtxV_neuroCHARGE/BrainPhenoNames.txt")
BrainPhenoNames <- read.csv(text = BrainPhenoNames, header=F)
head(BrainPhenoNames)

setdiff(names(d_brain_ICV),BrainPhenoNames$V1)
setdiff(BrainPhenoNames$V1,names(d_brain_ICV))
names(d_brain_ICV) = str_replace(names(d_brain_ICV),'totalsurface','WhiteSurfArea')
names(d_brain_ICV) = str_replace(names(d_brain_ICV),'thickness_globalmeanmean','MeanThickness_thickness')
names(d_brain_ICV) = str_replace(names(d_brain_ICV),"wholebrain_estimatedtotalintracranial_volume","ICV")

d_brain_ICV$lh_temporalpole_area <- NA
d_brain_ICV$lh_temporalpole_thickness <- NA
d_brain_ICV$lh_temporalpole_volume <- NA
d_brain_ICV$rh_temporalpole_area <- NA
d_brain_ICV$rh_temporalpole_thickness <- NA
d_brain_ICV$rh_temporalpole_volume <- NA

col_names = c('eid',BrainPhenoNames$V1)
d_brain_ICV = d_brain_ICV %>% dplyr::select(all_of(col_names))
identical(names(d_brain_ICV),col_names)
head(d_brain_ICV %>% dplyr::select(matches("temporalpole")))
# write_tsv(d_brain_ICV %>% dplyr::select(matches("eid|thickness")), 
#           "~/Downloads/d_brain_thickness.tsv")

d_T2D = fread("/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb/data/SummaryDiabetesDiagnosis_baseline.txt")

write_tsv(d_brain_ICV,"/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts/ukb677610_d_brain_ICV.txt")
write_tsv(d_geno,"/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts/ukb677610_d_geno.txt")
write_tsv(d_bmi_covs %>% left_join(d_T2D),
          "/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts/ukb677610_d_bmi_covs.txt")

# create an html RMD file by rendering ----------------------------------------
# rmarkdown::render("~/Documents/scripts/neuroCHARGE_scripts/Obesity_CtxV/format_ukb677610_data.R")
library(ggplot2)
df = d_brain_ICV %>% dplyr::select(matches("rh_")) %>% 
  dplyr::select(matches("thickness"))
names(df) = str_remove(names(df))
ggplot(stack(), aes(x = ind, y = values)) +
  geom_boxplot()

ggplot(stack(d_brain_ICV %>% dplyr::select(matches("lh_")) %>% 
               dplyr::select(matches("thickness"))), aes(x = ind, y = values)) +
  geom_boxplot()
