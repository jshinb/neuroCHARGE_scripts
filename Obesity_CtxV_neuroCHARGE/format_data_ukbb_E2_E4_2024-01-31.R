#' ---
#' title: "Formatting UK Biobank data"
#' author: "Jean Shin"
#' date: "2024-01-31"
#' 

#' Need to write a function to generate APOE4 status 
#' This VERSION re-defined infer_APOE() function
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

infer_APOE = function(data, rs429358_name, rs7412_name){
  data$APOE=NA
  
  # assuming C/T are minor/major allele for rs429358
  ind1 = !is.na(data$rs429358_chr19_C_T);print(sum(ind1))
  ind1_CC= ind1 & data$rs429358_chr19_C_T=="CC";print(sum(ind1_CC))#757
  ind1_CT= ind1 & data$rs429358_chr19_C_T=="CT";print(sum(ind1_CT))#8574
  ind1_TT= ind1 & data$rs429358_chr19_C_T=="TT";print(sum(ind1_TT))#23742
  
  # assuming T/C are minor/major allele for rs7412
  ind2 = !is.na(data$rs7412_chr19_T_C);print(sum(ind2))
  ind2_TT= ind2 & data$rs7412_chr19_T_C=="TT";print(sum(ind2_TT))#193
  ind2_TC= ind2 & data$rs7412_chr19_T_C=="TC";print(sum(ind2_TC))#4950
  ind2_CC= ind2 & data$rs7412_chr19_T_C=="CC";print(sum(ind2_CC))#27930
  
  data$APOE = NA
  #
  data$APOE[ind1_CC & ind2_TT] = 'E1E1';print(sum(ind1_CC & ind2_TT))#0
  data$APOE[ind1_CC & ind2_TC] = 'E1E4';print(sum(ind1_CC & ind2_TC))#0
  data$APOE[ind1_CC & ind2_CC] = 'E4E4';print(sum(ind1_CC & ind2_CC))#751
  #
  data$APOE[ind1_CT & ind2_TT] = 'E1E2';print(sum(ind1_CT & ind2_TT))#0
  data$APOE[ind1_CT & ind2_TC] = 'E2E4/E1E3';print(sum(ind1_CT & ind2_TC))#804
  data$APOE[ind1_CT & ind2_CC] = 'E3E4';print(sum(ind1_CT & ind2_CC))#7,683
  #
  data$APOE[ind1_TT & ind2_TT] = 'E2E2';print(sum(ind1_TT & ind2_TT))#189
  data$APOE[ind1_TT & ind2_TC] = 'E2E3';print(sum(ind1_TT & ind2_TC))#4,083
  data$APOE[ind1_TT & ind2_CC] = 'E3E3';print(sum(ind1_TT & ind2_CC))#19,212
  
  ## ambiguous haplotypes are excluded from the analyses 
  data$E4status <- NA
  data$E4status[!is.na(data$APOE) & data$APOE %in% c('E3E4', 'E4E4')] <- 'e4-carrier'
  data$E4status[!is.na(data$APOE) & data$APOE %in% c('E2E2', 'E2E3')] <- 'e4-noncarrier'

  data$APOE_simple = NA
  data$APOE_simple[!is.na(data$APOE) & data$APOE %in% c('E3E4', 'E4E4')] <- 'E4'
  data$APOE_simple[!is.na(data$APOE) & data$APOE %in% c('E2E2', 'E2E3')] <- 'E2'
  data$APOE_simple[!is.na(data$APOE) & data$APOE == 'E3E3'] <- 'E3'
  
  data
}

# 
rois = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/neuroCHARGE_WMH_ctxTH_age/neuroCHARGE_tests/freesurfer_34rois.txt',header=F)$V1
ctx_col_names = c(paste(rois,"right","hemisphere",sep="_"),
                  paste(rois,"left","hemisphere",sep="_"))

#Read in-house FreeSurfer-processed values ------------------------------------
d_freesurfer=fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb/data/freesurfer_combined.txt')
d_other_data = fread('ukb_data/ukbb41448_other_data_2023-01-11.tsv')
head(d_other_data)

d_freesurfer = data.frame(d_freesurfer)
select_cols = c(paste("rh",rois,"volume",sep="_"),
                paste("lh",rois,"volume",sep="_"),
                paste('rh',rois,'area',sep="_"),
                paste('lh',rois,'area',sep="_"),
                'rh_WhiteSurfArea_area',
                'lh_WhiteSurfArea_area',
                paste('rh',rois,'thickness',sep="_"),
                paste('lh',rois,'thickness',sep="_"),
                'rh_MeanThickness_thickness',
                'lh_MeanThickness_thickness')

d_merge = d_other_data %>% 
  right_join(d_freesurfer%>% 
               dplyr::select(all_of(c("ID",select_cols))) %>% dplyr::rename(eid=ID)) %>% 
  dplyr::rename(mri_site=uk_biobank_assessment_centre_f54_2_0) %>% 
  mutate(mri_site = factor(mri_site))
d_merge = d_merge %>% 
  right_join(subset(d_freesurfer,select=c(ID,aseg_volume_EstimatedTotalIntraCranialVol))%>% 
            dplyr::rename(eid = ID, ICV=aseg_volume_EstimatedTotalIntraCranialVol))
(dim(d_merge))

mri_site_coding = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb/data/ukb_MRI_centre_code_2020-11-30.tsv')
d_merge = d_merge %>% mutate(mri_site_coded = NA)
for(site_id in unique(d_merge$mri_site)){
  d_merge$mri_site_coded[d_merge$mri_site==site_id] = mri_site_coding$center[mri_site_coding$code==site_id]
}
d_merge = d_merge %>% mutate(mri_site_coded = str_split(d_merge$mri_site_coded," ",simplify = T)[,1])
d_merge = d_merge %>% dplyr::select(-mri_site) %>% dplyr::rename(mri_site=mri_site_coded)

# rm(d_ctx_area,d_ctx_thickness,d_ctx_volume,d_other_data);gc(T)

# other data2 -----------------------------------------------------------------
d_other_data2 = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/neuroCHARGE_obesity_ctxTH_metabolomics/data/ukbb41448_bmi_sbp_ukb41450_cholesterol_2023-09-26.tsv')
d_merge = d_merge %>% left_join(d_other_data2 %>% dplyr::select(eid,body_mass_index_bmi_f21001_0_0))

#exclusions -------------------------------------------------------------------
library(ukbtools)
## dropouts
DropOut_subs <- fread("~/OneDrive - SickKids/ukbb/data/w43688_2023-04-25.csv")$V1
d_merge = subset(d_merge,!eid %in% DropOut_subs)

## genetically validated british whites 
d_white = fread('~/OneDrive - SickKids/ukbb/data/ukbb_genetic_white.txt')
d_merge = d_merge %>% filter(eid %in% d_white$eid[!is.na(d_white$genetic_white)])#33109

## related 
kinship.info <- fread("/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb/data/ukb43688_rel_s488264.dat")
kinship.info <- subset(kinship.info, ID1>0)
kinship.info = subset(kinship.info,Kinship>0)
kinship.anal = subset(kinship.info,ID1 %in% d_merge$eid & kinship.info$ID2 %in% d_merge$eid)#452

remove.ids <- ukb_gen_samples_to_remove(kinship.anal, d_merge$eid, cutoff = 0.0884)#420#(default 0.0884 includes pairs with greater than 3rd-degree relatedness)
d_merge = d_merge %>% filter(!eid %in% remove.ids)
dim(d_merge)#33,125

d_genetic_sex = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb/data/ukbb41449_genetic_sex_2023-10-20.tsv')
d_genetic_sex = d_genetic_sex %>% dplyr::rename(genetic_sex=genetic_sex_f22001_0_0)

d_merge = d_merge %>% 
  dplyr::rename(age=age_when_attended_assessment_centre_f21003_2_0)%>%#age at MRI
  dplyr::rename(sex=sex_f31_0_0,
                weight=weight_f21002_2_0,
                height=height_f12144_2_0,
                bmi0=body_mass_index_bmi_f21001_0_0,
                bmi1=body_mass_index_bmi_f21001_2_0) %>%#BMI - measured at baseline
  mutate(age.group=cut(age,breaks=c(35,45,55,65,75,85))) %>%
  left_join(d_genetic_sex)
# use genetic sex instead (20 inconsistent values)
d_merge = d_merge %>% dplyr::rename(reported_sex = sex) %>% dplyr::rename(sex=genetic_sex)
get_edu = function(x){
  high.edu = c('A levels/AS levels or equivalent','College or University degree','Other professional qualifications eg: nursing, teaching')
  if(any(x %in% high.edu)){
    y = 'high'
  }else{
    if(any(na.omit(x) == 'Prefer not to answer')){
      y <- NA
    }else{
      y = 'low'
    }
  }
}
d_merge = d_merge %>%
  mutate(edu = apply(data.frame(d_merge[,paste("qualifications_f6138_0",0:5,sep="_")]),1,get_edu))
head(subset(d_merge,select=c(paste("qualifications_f6138_0",0:5,sep="_"),'edu'),edu=="low"))

table(d_merge$age.group)
# (35,45] (45,55] (55,65] (65,75] (75,85] 
#      2    5649   12546   13392    1536 

agecat2 = rownames(table(d_merge$age.group,d_merge$sex))
table(d_merge$sex,d_merge$age.group,useNA='a')
#        (35,45] (45,55] (55,65] (65,75] (75,85] <NA>
# Female       2    3198    7063    6466     628    0
# Male         0    2451    5483    6926     908    0
# <NA>         0       0       0       0       0    0
d_merge$age.group[d_merge$age.group=="(35,45]"] <- "(45,55]"#
names(d_merge)

# stroke and dementia ---------------------------------------------------------
d_stroke = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb_insulin_resistance/data/ukbbStroke_Dementia_at_anytime.txt')
dim(d_merge %>% filter(eid %in% d_stroke$eid))#352 individuals with stroke/dementia
d_merge = d_merge %>% filter(!eid %in% d_stroke$eid)#32773


#extracting brain_pheno -------------------------------------------------------
rois = c("bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", 
         "cuneus", "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal", 
         "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal", 
         "lingual", "medialorbitofrontal", "middletemporal", "parahippocampal", 
         "paracentral", "parsopercularis", "parsorbitalis", "parstriangularis", 
         "pericalcarine", "postcentral", "posteriorcingulate", "precentral", 
         "precuneus", "rostralanteriorcingulate", "rostralmiddlefrontal", 
         "superiorfrontal", "superiorparietal", "superiortemporal", "supramarginal", 
         "frontalpole", "temporalpole", "transversetemporal", "insula")
varnames = 'eid'
varnames = c(varnames,paste("rh",rois,"volume",sep="_"),
             paste("lh",rois,"volume",sep="_"),
             paste('rh',rois,'area',sep="_"),
             paste('lh',rois,'area',sep="_"),
             'rh_WhiteSurfArea_area',
             'lh_WhiteSurfArea_area',
             paste('rh',rois,'thickness',sep="_"),
             paste('lh',rois,'thickness',sep="_"),
             'rh_MeanThickness_thickness',
             'lh_MeanThickness_thickness')
varnames = c(varnames,'ICV')

d_brain = d_merge %>% dplyr::select(all_of(varnames))

# extracting genotype and genetic PC data -------------------------------------
d_geno = d_merge %>% dplyr::select(eid,contains("genoPC"),contains("APOE"))
head(d_geno)
table(d_geno$APOE)

# extracting covariates -------------------------------------------------------
varnames_bmi_covariates = c('eid',
                            'bmi1',# bmi@ MRI scan
                            'age',
                            'age_when_attended_assessment_centre_f21003_0_0',
                            'sex',
                            'edu',
                            'mri_site')
d_bmi_covariates = d_merge %>% dplyr::select(all_of(varnames_bmi_covariates)) %>%
  dplyr::rename(bmi = bmi1) %>%
  dplyr::rename(age_bmi = age_when_attended_assessment_centre_f21003_0_0) %>%
  dplyr::rename(age_mri = age)
head(d_bmi_covariates)

# 
# 1. write 3 files ---------------------------------------------------------------
write_tsv(d_brain,"/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts/BrainData.txt")
write_tsv(d_geno,"/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts/GenoData.txt")#APOE + genoPC
write_tsv(d_bmi_covariates %>% dplyr::select(eid,mri_site),"/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts/CovData_UKBB.txt")#APOE + genoPC

# 
d_merge = d_brain %>%   
  left_join(d_bmi_covariates) %>%
  left_join(d_geno)
head(d_merge)
names(d_merge)

# additional covairates: smoking, bp, medication  -----------------------------
d_additional_covs = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb/data/ukbb41448_41450_sbp_dbp_smoking_medication_covariates.tsv')
names(d_additional_covs)
d_additional_covs_formatted = d_additional_covs %>% dplyr::select(eid)

##SBP
SBP = d_additional_covs %>% 
  dplyr::select(eid,
                systolic_blood_pressure_automated_reading_f4080_2_0,
                systolic_blood_pressure_automated_reading_f4080_2_1,
                systolic_blood_pressure_manual_reading_f93_2_0,
                systolic_blood_pressure_manual_reading_f93_2_1)
head(SBP)
names(SBP)[-1] <- c("auto0","auto1","manu0","manu1")
SBP = SBP %>% mutate(sbp = apply(data.frame(subset(SBP,select=-eid)),1,mean,na.rm=T))
SBP %>% filter(!is.na(sbp))
cmatSBP = cor(SBP %>% dplyr::select(-eid),use='p')

##DBP
DBP = d_additional_covs %>% 
  dplyr::select(eid,
                diastolic_blood_pressure_automated_reading_f4079_2_0,
                diastolic_blood_pressure_automated_reading_f4079_2_1,
                diastolic_blood_pressure_manual_reading_f94_2_0,
                diastolic_blood_pressure_manual_reading_f94_2_1)
head(DBP)
names(DBP)[-1] <- c("auto0","auto1","manu0","manu1")
DBP = DBP %>% mutate(dbp = apply(data.frame(subset(DBP,select=-eid)),1,mean,na.rm=T))
DBP %>% filter(!is.na(dbp))
cmatDBP = cor(DBP %>% dplyr::select(-eid),use='p')

## HTN medication - at baseline
HTNmed = d_additional_covs %>%
  dplyr::select(eid,
                medication_for_cholesterol_blood_pressure_or_diabetes_f6177_2_0,
                medication_for_cholesterol_blood_pressure_or_diabetes_f6177_2_1,
                medication_for_cholesterol_blood_pressure_or_diabetes_f6177_2_2,#male
                medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_f6153_2_0,
                medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_f6153_2_1,
                medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_f6153_2_2,
                medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_f6153_2_3)

names(HTNmed)[-1] = c("medM0","medM1","medM2","medF0","medF1","medF2","medF3")

take.htn.med = function(x){
  names(x) = NULL
  if(any(!is.na(x)) & any(x %in% "Blood pressure medication")){
    y=1
    } else if(any(!is.na(x)) & !any(x %in% "Blood pressure medication")){
      y=0
    }else{
      y=NA
    }
  y
}
HTNmed = HTNmed %>% mutate(htnmed = apply(subset(HTNmed,select = -eid),1,take.htn.med))

## derive hypertension status
HTNstatus = SBP %>% dplyr::select(eid,sbp) %>%
  left_join(DBP %>% dplyr::select(eid,dbp)) %>%
  left_join(HTNmed %>% dplyr::select(eid,htnmed))

HTNstatus = HTNstatus %>% mutate(HTNstatus = ifelse((sbp>=140|dbp>=90)|htnmed==1,1,0))
ind = !is.na(HTNstatus$sbp) & !is.na(HTNstatus$dbp);print(sum(ind))
ind = ind & (HTNstatus$sbp<140&HTNstatus$dbp<90);print(sum(ind))
ind = ind & is.na(HTNstatus$htnmed);print(sum(ind))
ind = ind & is.na(HTNstatus$htnmed);print(sum(ind))
  
HTNstatus %>% filter(ind)
HTNstatus$HTNstatus[ind] = 0
table(HTNstatus$HTNstatus)
HTNstatus = HTNstatus %>% filter(!is.na(HTNstatus))

## smoking -yes/no or current/former/never
smoking = d_additional_covs %>% 
  dplyr::select(eid,smoking_status_f20116_2_0)
names(smoking)[2] = 'current_smoking'
smoking %>% filter(!is.na(current_smoking))
table(smoking$current_smoking)

## merge to the original cov data
d_bmi_covariates = d_bmi_covariates %>% 
  left_join(HTNstatus, by='eid') %>% 
  left_join(smoking)

# type 2 diabetes -------------------------------------------------------------
d_T2D = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb/data/SummaryDiabetesDiagnosis_baseline.txt')

## merge to the original cov data
d_bmi_covariates = d_bmi_covariates%>%
  left_join(d_T2D)

# smoking status - redefinition -----------------------------------------------
table(d_bmi_covariates$current_smoking)
current_smoking_yes = d_bmi_covariates[['current_smoking']]
current_smoking_yes[!is.na(current_smoking_yes) & current_smoking_yes=="Prefer not to answer"] <- NA

ind.not.na = !is.na(current_smoking_yes);print(sum(ind.not.na))
ind = ind.not.na & current_smoking_yes == "Current";print(sum(ind))
ind.not.current = ind.not.na & current_smoking_yes != "Current";print(sum(ind.not.current))
current_smoking_yes[ind] = 'Yes'
current_smoking_yes[ind.not.current] = 'No'
d_bmi_covariates = d_bmi_covariates %>% mutate(current_smoking_status = current_smoking_yes)
head(d_bmi_covariates)

# write out the COV and APOE4 geno files --------------------------------------
data_dir = '/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts'
psych::describe(d_bmi_covariates)
write_tsv(d_bmi_covariates %>% dplyr::select(-mri_site),
          file.path(data_dir,"BMI_and_CovData.txt"))#APOE + genoPC
# genotype PC data ------------------------------------------------------------ 
d_genoPC=fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb_insulin_resistance/data/ukbb41449_genoPC_2023-04-10.tsv')
head(d_genoPC)

# APOE genotypes --------------------------------------------------------------
d_APOE = fread('ukb_data/hardcalled_genotypes_for_APOE_genotypes_2023-05-04.tsv')
d_APOE = d_APOE %>% filter(IID %in% d_merge$eid)#32722

head(d_APOE)
infer_APOE = function(data, rs429358_name, rs7412_name){
  data$APOE=NA
  
  # assuming C/T are minor/major allele for rs429358
  ind1 = !is.na(data$rs429358_chr19_C_T);print(sum(ind1))
  ind1_CC= ind1 & data$rs429358_chr19_C_T=="CC";print(sum(ind1_CC))#757
  ind1_CT= ind1 & data$rs429358_chr19_C_T=="CT";print(sum(ind1_CT))#8574
  ind1_TT= ind1 & data$rs429358_chr19_C_T=="TT";print(sum(ind1_TT))#23742
  
  # assuming T/C are minor/major allele for rs7412
  ind2 = !is.na(data$rs7412_chr19_T_C);print(sum(ind2))
  ind2_TT= ind2 & data$rs7412_chr19_T_C=="TT";print(sum(ind2_TT))#193
  ind2_TC= ind2 & data$rs7412_chr19_T_C=="TC";print(sum(ind2_TC))#4950
  ind2_CC= ind2 & data$rs7412_chr19_T_C=="CC";print(sum(ind2_CC))#27930
  
  data$APOE = NA
  #
  data$APOE[ind1_CC & ind2_TT] = 'E1E1';print(sum(ind1_CC & ind2_TT))#0
  data$APOE[ind1_CC & ind2_TC] = 'E1E4';print(sum(ind1_CC & ind2_TC))#0
  data$APOE[ind1_CC & ind2_CC] = 'E4E4';print(sum(ind1_CC & ind2_CC))#751
  #
  data$APOE[ind1_CT & ind2_TT] = 'E1E2';print(sum(ind1_CT & ind2_TT))#0
  data$APOE[ind1_CT & ind2_TC] = 'E2E4/E1E3';print(sum(ind1_CT & ind2_TC))#804
  data$APOE[ind1_CT & ind2_CC] = 'E3E4';print(sum(ind1_CT & ind2_CC))#7,683
  #
  data$APOE[ind1_TT & ind2_TT] = 'E2E2';print(sum(ind1_TT & ind2_TT))#189
  data$APOE[ind1_TT & ind2_TC] = 'E2E3';print(sum(ind1_TT & ind2_TC))#4,083
  data$APOE[ind1_TT & ind2_CC] = 'E3E3';print(sum(ind1_TT & ind2_CC))#19,212
  
  ## ambiguous haplotypes are excluded from the analyses 
  #data$E4status <- NA
  #data$E4status[!is.na(data$APOE) & data$APOE %in% c('E3E4', 'E4E4')] <- 'e4-carrier'
  #data$E4status[!is.na(data$APOE) & data$APOE %in% c('E2E2', 'E2E3','E3E3')] <- 'e4-noncarrier'

  data$E4status <- NA
  data$E4status[!is.na(data$APOE) & data$APOE %in% c('E3E4', 'E4E4')] <- 'e4-carrier'
  data$E4status[!is.na(data$APOE) & data$APOE %in% c('E2E2', 'E2E3')] <- 'e4-noncarrier' #no e3-e3

  data$APOE_simple = NA
  data$APOE_simple[!is.na(data$APOE) & data$APOE %in% c('E3E4', 'E4E4')] <- 'E4'
  data$APOE_simple[!is.na(data$APOE) & data$APOE %in% c('E2E2', 'E2E3')] <- 'E2'
  data$APOE_simple[!is.na(data$APOE) & data$APOE == 'E3E3'] <- 'E3'
  
  data
}

d_APOE_wiE4 = infer_APOE(d_APOE, rs429358_name="rs429358_chr19_C_T", rs7412_name="rs7412_chr19_T_C")

write_tsv(d_APOE_wiE4 %>% dplyr::rename(eid=IID) %>% left_join(d_genoPC,by="eid"),file.path(data_dir,"APOE_and_genoPCs_E2_vs_E4.txt"))

# create an htmp RMD file by rendering ----------------------------------------
# rmarkdown::render("~/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts/format_data_ukbb.R")
setwd('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts')
