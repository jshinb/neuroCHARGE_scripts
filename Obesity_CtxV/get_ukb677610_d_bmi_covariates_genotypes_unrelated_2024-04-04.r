# Select individuals ----------------------------------------------------------
# unrelated
# White British
# no history of stroke and dementia
# 
#  module load gcc/8.3.0 intel/2019u4 r/4.1.2;R

# load libraries and functions ------------------------------------------------
library(ukbtools)
library(stringr)
library(dplyr)
library(data.table)
library(readr)

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
  data$E4status[!is.na(data$APOE) & data$APOE %in% c('E2E2', 'E2E3','E3E3')] <- 'e4-noncarrier'
  data$APOE_simple = NA
  data$APOE_simple[!is.na(data$APOE) & data$APOE %in% c('E3E4', 'E4E4')] <- 'E4'
  data$APOE_simple[!is.na(data$APOE) & data$APOE %in% c('E2E2', 'E2E3')] <- 'E2'
  data$APOE_simple[!is.na(data$APOE) & data$APOE == 'E3E3'] <- 'E3'
  
  data
}

setwd("/scratch/t/tpaus/jshinb/KEEP/UKBB/datasets/ukb677610/")

# read in data ------------------------------------------------------------------------
d = ukb_df('ukb677610')

## other covariates
# bmi - f21001f21003
# age at recruitment - 21022_0_0 
# age when attended assessment centre - 21003_X_X
# sex = genetic sex 
# education 
# sbp - 93-0.0, 93-0.1 (average of two measurements) - auto
# dbp - 94-0.0, 94-0.1 (average of two measurements) - auto
# sbp_auto - 4080-0.0, 4080-0.1
# dbp_auto - 4079-0.0, 4089-0.1
# htn medication - 6153_0_1, 6153_0_2, 6153_0_3 [coding:2, NA: -7, -1, -3]
# htn status - based on SBP, DBP and HTN medication 
# current_smoking - 20116-0.0
# from '~/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts/format_data_ukbb.R'
d_brain =fread('/scratch/t/tpaus/jshinb/KEEP/UKBB_JS/ukb677610_d_brain_unrelated.txt')#46,193

d_additional_covs_formatted = d %>% dplyr::select(eid)

## blood pressure and hypertension status -------------------------------------
SBP = d %>% filter(eid %in% d_brain$eid) %>%
  dplyr::select(eid,
                systolic_blood_pressure_automated_reading_f4080_2_0,
                systolic_blood_pressure_automated_reading_f4080_2_1,
                systolic_blood_pressure_manual_reading_f93_2_0,
                systolic_blood_pressure_manual_reading_f93_2_1)
head(SBP)
names(SBP)[-1] <- c("auto0","auto1","manu0","manu1")
SBP = SBP %>% mutate(sbp = apply(data.frame(subset(SBP,select=-eid)),1,mean,na.rm=T))
print(dim(SBP))
cmatSBP = cor(SBP %>% dplyr::select(-eid),use='p')

##DBP
DBP = 
  d %>% filter(eid %in% d_brain$eid) %>%
  dplyr::select(eid,
                diastolic_blood_pressure_automated_reading_f4079_2_0,
                diastolic_blood_pressure_automated_reading_f4079_2_1,
                diastolic_blood_pressure_manual_reading_f94_2_0,
                diastolic_blood_pressure_manual_reading_f94_2_1)
head(DBP)
names(DBP)[-1] <- c("auto0","auto1","manu0","manu1")
DBP = DBP %>% mutate(dbp = apply(data.frame(subset(DBP,select=-eid)),1,mean,na.rm=T))
cmatDBP = cor(DBP %>% dplyr::select(-eid),use='p')

## HTN medication - at baseline
HTNmed = 
d %>% filter(eid %in% d_brain$eid) %>%
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
ind = !is.na(HTNstatus$sbp) & !is.na(HTNstatus$dbp);print(sum(ind))#39762
ind = ind & (HTNstatus$sbp<140&HTNstatus$dbp<90);print(sum(ind))#20882
# ind = ind & is.na(HTNstatus$htnmed);print(sum(ind))
# HTNstatus$HTNstatus[ind] = 0 # don't know whether they are on HTN medication and control BP
table(HTNstatus$HTNstatus,useNA='a')#5225 with no medication information


## smoking -yes/no or current/former/never ------------------------------------
smoking = 
  d %>% filter(eid %in% d_brain$eid) %>%
  dplyr::select(eid,smoking_status_f20116_2_0)
names(smoking)[2] = 'current_smoking'
smoking %>% filter(!is.na(current_smoking))
table(smoking$current_smoking)

current_smoking_yes = smoking[['current_smoking']]
print(table(current_smoking_yes,useNA='a'))

current_smoking_yes[!is.na(current_smoking_yes) & current_smoking_yes=="-3"] <- NA#Prefer not to answer
print(table(current_smoking_yes,useNA='a'))

current_smoking_yes = ifelse(current_smoking_yes==2, "Yes", "No")
print(table(current_smoking_yes,useNA='a'))

current_smoking_yes
smoking = smoking %>% mutate(current_smoking = current_smoking_yes)
#    No   Yes  <NA>
# 42924  1499   467

## age, sex and bmi -----------------------------------------------------------
d_age_sex_bmi = 
  d %>% filter(eid %in% d_brain$eid) %>%
  dplyr::select(matches("eid|f21003"))#age
names(d_age_sex_bmi)[-1] = c('age0','age1','age2','age3')

d_age_sex_bmi = 
  d_age_sex_bmi %>% dplyr::select(-age1,-age3) %>%
  left_join( 
    d %>% 
      select(all_of(c("eid","sex_f31_0_0"))) %>%
      mutate(sex = ifelse(sex_f31_0_0==0, 'Female', "Male"))#0 female
  ) %>% 
  dplyr::select(-sex_f31_0_0)

d_age_sex_bmi = 
  d_age_sex_bmi %>% 
  left_join(
    d %>% 
      filter(eid %in% d_brain$eid) %>%
      dplyr::select(matches("eid|f21001"))#bmi
  ) %>%
  dplyr::select(-body_mass_index_bmi_f21001_1_0,-body_mass_index_bmi_f21001_3_0) %>%
  dplyr::rename(bmi0 = body_mass_index_bmi_f21001_0_0,bmi2 = body_mass_index_bmi_f21001_2_0) 

print(dim(d_age_sex_bmi))
print(head(d_age_sex_bmi))

## education ------------------------------------------------------------------
translate_edu_code = function(x){
  y = rep(NA,length(x))
  ind.notNA = !is.na(x)
  codings = c(-7, -3,  1,  2,  3,  4,  5,  6)
  meanings = 'None of the above
Prefer not to answer
College or University degree
A levels/AS levels or equivalent
O levels/GCSEs or equivalent
CSEs or equivalent
NVQ or HND or HNC or equivalent
Other professional qualifications eg: nursing, teaching'
  meanings = unlist(stringr::str_split(meanings,"\n"))
  i=0
  for (ci in codings){
    i=i+1
    y[ind.notNA & x==as.character(ci)] <- meanings[i]   
  }
  y
}
get_edu = function(x){
  #https://academic.oup.com/ije/article/51/3/885/6521336
  #Carter, Alice R et al. “Educational attainment as a modifier for the effect of polygenic scores for cardiovascular risk factors: cross-sectional and prospective analysis of UK Biobank.” International journal of epidemiology vol. 51,3 (2022): 885-897. doi:10.1093/ije/dyac002
  high.edu = c('College or University degree',
               'NVQ or HND or HNC or equivalent',
               'Other professional qualifications eg: nursing, teaching')
  
  low.edu = c('A levels/AS levels or equivalent',
              'O levels/GCSEs or equivalent',
              'CSEs or equivalent',
              'None of the above'
  )
  na.edu = c('Prefer not to answer')
  
 y = ifelse(x %in% high.edu, 'high', x)
 y = ifelse(x %in% low.edu, 'low', y)
 y = ifelse(x %in% na.edu, NA, y)
 y
}
get_edu_high_low_overall = function(x){
  if(any(sum(!is.na(x)))){
    y.high = any(x[!is.na(x)]=="high")
    y.low = all(x[!is.na(x)]=="low")
    y.NA <- FALSE
  }else{
    y.high <- FALSE
    y.low <- FALSE
    y.NA <- TRUE
  }
  y = c('high','low',NA)[c(y.high,y.low,y.NA)]
  y
}

d_education = d %>% filter(eid %in% d_brain$eid) %>%
  dplyr::select(matches("eid|f6138_"))
d_edu_coding = fread('/scratch/t/tpaus/jshinb/KEEP/UKBB_JS/education_coding.txt')
head(d_edu_coding)

tmp_edu = data.frame(d_education %>% dplyr::select(-eid)) 
tmp_edu = apply(tmp_edu,2,translate_edu_code)
head(tmp_edu)

edu_high_low = c()
for( i in 1:ncol(tmp_edu)){
  edu_high_low =cbind(edu_high_low, get_edu(tmp_edu[,i]))
}
colnames(edu_high_low) <- colnames(tmp_edu) <- 
  str_replace(colnames(tmp_edu),"qualifications_f6138_","edu_")

edu = data.frame(edu_high_low) %>% dplyr::select(matches("edu_0_|edu_2_"))
edu_0 = apply(edu[,1:6],1,get_edu_high_low_overall)
edu_2 = apply(edu[,7:12],1,get_edu_high_low_overall)
table(edu_2, edu_0,useNA='a')

# some of them had 'high' education at t0, but lower education level at t2, vice versa
# I will assume t1 information
ind = edu_0 == 'high' & edu_2 == "low"
ind = !is.na(ind) & ind; print(sum(ind))
head(tmp_edu[ind])

d_education = d_education %>% dplyr::select(eid) %>% mutate(edu0 = edu_0, edu2=edu_2)

# merge_covariates ------------------------------------------------------------
d_bmi_covariates = d_age_sex_bmi %>%
  left_join(HTNstatus) %>%
  left_join(smoking) %>%
  left_join(d_education)
readr::write_tsv(d_bmi_covariates,
                 "/scratch/t/tpaus/jshinb/KEEP/UKBB_JS/ukb677610_d_bmi_covariates.txt")
# genetic information ---------------------------------------------------------
d_apoE = fread('/scratch/t/tpaus/jshinb/KEEP/UKBB_JS/hardcalled_genotypes_for_APOE_genotypes_2023-05-04.tsv') %>%
  filter(IID %in% d_brain$eid) %>% #487409 - 44890
  dplyr::rename(eid=IID)#44,890

head(d_apoE)
d_apoE4 = infer_APOE(d_apoE)
round(prop.table(table(d_apoE4$APOE))*100,3)
#E1E4      E2E2      E2E3 E2E4/E1E3      E3E3      E3E4      E4E4
#   1       264      5582      1051     26579     10421       992

# E1E4      E2E2      E2E3 E2E4/E1E3      E3E3      E3E4      E4E4
# 0.002     0.588    12.435     2.341    59.209    23.215     2.210

# from:https://www.nature.com/articles/s41398-024-02848-5/tables/1
round(c("E2-carrier"=3842,"E3E3"=17206,"E3E4"=6794,"E4E4"=652)/sum(c("E2-carrier"=3842,"E3E3"=17206,"E3E4"=6794,"E4E4"=652))*100,3)
#E2-carrier       E3E3       E3E4       E4E4 
#.   13.484     60.385     23.844      2.288 

d_genoPC = d %>% dplyr::select(matches("eid|f22009"))
names(d_genoPC) = str_replace(names(d_genoPC),
                              "genetic_principal_components_f22009_0_","genoPC")
d_genoPC = d_genoPC %>% filter(eid %in% d_brain$eid)
d_apoE4_genoPC = d_apoE4 %>% left_join(d_genoPC)

readr::write_tsv(d_apoE4_genoPC,"/scratch/t/tpaus/jshinb/KEEP/UKBB_JS/ukb677610_d_apoE4_genoPC.txt")
