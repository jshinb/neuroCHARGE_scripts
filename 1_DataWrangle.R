#*****************************************************************************#
#
# Step 1: Wrangle data 
#
#*****************************************************************************#
cat("\nStep 1: Starting reading and merging data files\n")

# standardized column names
id_name = IID
brainpheno_names = fread("BrainOutcome_Headers.txt",header = F)$V1#fixed
covariate_names = c(age,sex,current_smoking,hypertension,type2diabetes,ICVorBrainVolume,
                    setdiff(names(fread(CohortSpecificCovariates_File)),IID))

# d will include all the individuals with complete information on covariates
d = fread(BrainOutcome_File)
# Checking brain outcome column names ----------------------------------------
cat("Checking brain outcome data column names \n") 
isOK.brain_columns = (setdiff(names(d),brainpheno_names) == IID & all(brainpheno_names %in% names(d)))
if(!all(isOK.brain_columns)) {# not OK
  msg = cat("The following column name(s) should be corrected:\n",
            "-------------------------------------------------",
            setdiff(names(d),brainpheno_names)[setdiff(names(d),brainpheno_names)!=IID],
            "-------------------------------------------------",
            "\n**Please make sure the column names of your brain outcomes in \'BrainOutcome_File\' are the same as those listed in \'BrainOutcome_Headers.txt\'.**\n",
            sep="\n")
  stop(msg)
}

d = merge(d,na.omit(fread(Covariates_File)),by=IID)
if(!is.null(CohortSpecificCovariates_File)){
  d = merge(d,na.omit(fread(CohortSpecificCovariates_File)),by=IID)
}

# Checking and formatting covariate column names ----------------------------------------
cat("Checking and formatting column names \n") 

# name vectors
names(id_name) = "IID"
names(brainpheno_names) = brainpheno_names
names(covariate_names) = c("age","sex","current_smoking","hypertension","type2diabetes","ICVorBrainVolume",#required cov
                           setdiff(names(fread(CohortSpecificCovariates_File)),IID))# study-specific
sel_colnames = c(id_name,brainpheno_names,covariate_names)
# change user-specified column names to the 'standard' names


if(!all(sel_colnames %in% names(d))){
  cat("----------------------------------\n")
  cat(setdiff(sel_colnames,names(d)),sep="\n")
  cat("----------------------------------\n")
  msg="Above column names are incorrectly specified.\nPlease correct your column name(s)."
  stop(msg)
}
d = subset(d,select=sel_colnames); names(d) = names(sel_colnames)
d = d %>% mutate(logWMH = log(1+WMH))
sex_tmp = ifelse(d$sex==code_female,"F","M")
current_smoking_tmp = ifelse(d$current_smoking == code_current_smoking_yes,"yes","no")
hypertension_tmp = ifelse(d$hypertension == code_hypertensive_yes,"yes","no")
type2diabetes_tmp = ifelse(d$type2diabetes == code_t2d_yes,"yes","no")

d$sex = factor(sex_tmp,levels=c("F","M"))
d$current_smoking = factor(current_smoking_tmp, levels=c("no","yes"))
d$hypertension = factor(hypertension_tmp, levels=c("no","yes"))
d$type2diabetes = factor(type2diabetes_tmp, levels=c('no','yes'))
rm(list=ls(pattern = '_tmp'))

cat("\nStep 1: Finishing reading and merging data files\n")
