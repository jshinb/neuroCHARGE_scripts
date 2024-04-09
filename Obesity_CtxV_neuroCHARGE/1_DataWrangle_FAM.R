#*****************************************************************************#
#
# 1: Wrangle data 
#
#*****************************************************************************#
cat("\n1. Starting reading and merging data files\n")

# start here ------------------------------------------------------------------
## column names
id_name = IID

brainpheno_names = fread("BrainPhenoNames.txt",header = F)$V1#need to fix this

cohort_specific_cov_names = c()
if(!is.na(CohortSpecificCovariates_File)){
  cohort_specific_cov_names = names(fread(CohortSpecificCovariates_File))
}

# read data files -------------------------------------------------------------
## d will include all the individuals with complete information on covariates
d_brain = fread(BrainPheno_File)

## Checking column names of the brain data file -------------------------------
cat("Checking brain outcome data column names \n")
## brain data file must consist of 210 columns: 
## 1 (ID) + 68 (volume) + 70 (surface area) + 70 (thickness) + 1 (ICV)
diff_from_std_brain_names = setdiff(names(d_brain),brainpheno_names)
not_in_data = setdiff(brainpheno_names,names(d_brain))
if(!all(diff_from_std_brain_names==IID)) {# not OK **TO-DO** NEED TO CHECK
  msg = cat("The following column name(s) should be corrected:\n",
            "-------------------------------------------------",
            setdiff(names(d_brain),brainpheno_names)[setdiff(names(d_brain),brainpheno_names)!=IID],
            "-------------------------------------------------",
            "\n**Please make sure the column names of your brain outcomes in \'BrainPheno_File\' are the same as those listed in \'BrainPhenoNames.txt\'.**\n",
            sep="\n")
  stop(msg)
}
if(length(not_in_data)>0) {# not OK **TO-DO** NEED TO CHECK
  msg = cat("The following brain phenotypes are missing:\n",
            "-------------------------------------------------",
            not_in_data,
            "-------------------------------------------------",
            "\n**Please add the missing brain outcomes to \'BrainPheno_File\'.**\n",
            sep="\n")
  stop(msg)
}

## calculate total volumes, areas, and mean thicknesses in/across the 34 regions 
d_totalVol = get_VolSums(d_brain,IID)
d_totalSa = get_AreaSums(d_brain,IID)
d_meanTh = get_ThicknessMeans(d_brain,IID)
brainpheno_names_total_mean = c(names(d_totalVol)[names(d_totalVol)!=IID],
                                names(d_totalSa)[names(d_totalSa)!=IID],
                                names(d_meanTh)[names(d_meanTh)!=IID])

## merge brain data with bmi and covariate data 
d = d_totalVol %>% left_join(d_totalSa) %>% left_join(d_meanTh)#brain-names already formatted

### merge bmi and covariates
covariate_names = c(#user-defined covariate names
  IID,
  FID,
  bmi, #exposure
  age_MRI,
  age_BMI,
  sex,
  education,
  #ICV, brain
  current_smoking,
  hypertension,
  type2diabetes#,
#  setdiff(cohort_specific_cov_names,IID)
)
d = d %>% 
  left_join (fread(BMI_Covariates_File) %>% dplyr::select(all_of(covariate_names)),by=IID) %>% 
  left_join(d_brain %>% dplyr::select(all_of(c(IID,"ICV"))),by=IID)

#### if ages at MRI and at BMI differ, add 'time' between MRI and BMI as a covariate
Time_bw_MRI_BMI = NULL
if(age_MRI != age_BMI){
  d[['time']] = d[[age_MRI]] - d[[age_BMI]]
  if(var(d$time,na.rm=T)>0){
    Time_bw_MRI_BMI = 'time'
  }
}else{
  d[["age_BMI"]] = d[[age_MRI]]
  age_BMI = "age_BMI"
}
covariate_names = c(covariate_names,Time_bw_MRI_BMI)

### merge cohort-specific covariates
if(!is.na(CohortSpecificCovariates_File)){
  cat('\nMerging with cohort-specific covariate data.\n')
  d = d %>% left_join(na.omit(fread(CohortSpecificCovariates_File)),by=IID)
}

## **WORK** merge genotypes and genotype PCs -------------------------------------------
d_geno =fread(APOE4_genoPC_File)
d_geno = subset(d_geno,select=c(IID,APOE4_status,genoPCs))

d = d %>% left_join(d_geno)

###renaming bmi and covariate column names for standardization
###will get errors if any column is absent from 
std_covariate_names = c(
  "BMI", #exposure
  "age",
  "age_BMI",
  "sex",
  "education",
  "ICV",
  "current_smoking",
  "hypertension",
  "type2diabetes",
  "E4_status",
  Time_bw_MRI_BMI,
  setdiff(cohort_specific_cov_names,IID)
)

d = d %>%  
  dplyr::rename(
    "IID"=IID,
    "FID"=FID,
    "BMI"=bmi,
    "age"=age_MRI,
    "age_BMI"=age_BMI,
    "sex"=sex,
    "education"=education,
    "ICV"=ICV,
    "current_smoking" = current_smoking,
    "hypertension"=hypertension,
    "type2diabetes"=type2diabetes,
    "E4_status"=APOE4_status)

## Checking and formatting covariate column names -----------------------------
cat("Checking and formatting column names \n") 
if(!all(std_covariate_names %in% names(d))){
  cat("----------------------------------\n")
  cat(setdiff(std_covariate_names,names(d)),sep="\n")
  cat("----------------------------------\n")
  msg="Above column names are incorrectly specified.\nPlease correct your column name(s)."
  stop(msg)
}
sex_tmp = ifelse(d$sex==code_female,"F","M")
current_smoking_tmp = ifelse(d$current_smoking == code_current_smoking_yes,"yes","no")
hypertension_tmp = ifelse(d$hypertension == code_hypertensive_yes,"yes","no")
type2diabetes_tmp = ifelse(d$type2diabetes == code_type2diabetes_yes,"yes","no")
e4status_tmp = ifelse(d$E4_status == code_e4_carrier,"E4-carrier","E4-noncarrier")

d = d %>% mutate(sex = factor(sex_tmp,levels=c("F","M")),
                 current_smoking = factor(current_smoking_tmp, levels=c("no","yes")),
                 hypertension = factor(hypertension_tmp, levels=c("no","yes")),
                 type2diabetes = factor(type2diabetes_tmp, levels=c('no','yes')),
                 E4_status = factor(e4status_tmp, levels=c('E4-noncarrier','E4-carrier')))
rm(list=ls(pattern = '_tmp'))

# age category ----------------------------------------------------------------
d = d %>% mutate(age_group = ifelse(age>=65,"age_gt_65","age_lt_65")) %>%
  mutate(age_group = factor(age_group,levels=c("age_lt_65","age_gt_65")))

# BMI category ----------------------------------------------------------------
d = d %>% mutate(BMI_group = ifelse(BMI>=25,"BMI_gt_25","BMI_lt_25")) %>%
  mutate(BMI_group = factor(BMI_group,levels=c("BMI_lt_25","BMI_gt_25")))

cat("\nFinishing reading and merging data files\n")
