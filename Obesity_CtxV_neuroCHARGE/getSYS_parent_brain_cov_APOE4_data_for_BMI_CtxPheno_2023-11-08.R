
# prep: load packages, define functions and set working directory -------------
if(!is.element("pacman",installed.packages()[,1])){
  install.packages("pacman")
}

pacman::p_load(tidyverse, data.table, psych, 
               patchwork, GGally,corrplot,hrbrthemes,
               tableone,FactoMineR,factoextra,
               ppcor, lsmeans, multcomp, ModelMetrics,
               #caret, 
               gridExtra, Hmisc, pastecs, testit)
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
  ret = list(counts.NA = counts.NA, cleaned_data=dat)
  ret
}

wd = '/Users/jshin/Library/CloudStorage/OneDrive-SickKids/neuroCHARGE_WMH_ctxTH_age/WMH_ctxTH/Example_SYS'
setwd(wd)
sys_data_dir ='/Users/jshin/Library/CloudStorage/OneDrive-SickKids/SYS_Database'

# read in data files ---------------------------------------------------------- 
f=file.path(sys_data_dir,'jean_shin_dataanalysis_completespssubjects_20201120_004249.txt')
d = fread(f)
head(d)
ind0=str_detect(names(d),'SPS-MRI-Freesurfer5.3.0-MRI-Aug2018');print(sum(ind0))
ind1=ind0 & str_detect(names(d),'thickness');print(sum(ind1))
ind1=ind1 | str_detect(names(d),'area');print(sum(ind1))
ind1=ind1 | str_detect(names(d),'volume');print(sum(ind1))
pheno.names = names(d)[ind1]#QC not relevant?

d_brain = subset(d,select=c('uniqueID',pheno.names))
names(d_brain) = str_remove(names(d_brain),'SPS-MRI-Freesurfer5.3.0-MRI-Aug2018.')
d_brain = data.frame(d_brain)
d_brain[,-1] = sapply(d_brain[,-1],as.numeric)
psych::describe(d_brain)[['n']]#602 without NA

rois = unique(str_split(names(d_brain),"_",simplify = T)[,2])
r1 = which(rois=="bankssts")
r34 = which(rois=='insula')
rois = rois[r1:r34]

extract_columns = c()
for(pheno in c("volume",'area','thickness')){
	for(hemi in c('lh','rh')){
		extract_columns = c(extract_columns,paste(hemi,rois,pheno,sep='_'))
	}#for(hemi)
}#for(pheno)

length(intersect(extract_columns,names(d_brain)))#all

extract_columns = c('uniqueID',extract_columns,
	'derived_lefthemvolume','derived_righthemvolume',
	'lh_WhiteSurfArea_area','rh_WhiteSurfArea_area',
	'lh_MeanThickness_thickness','rh_MeanThickness_thickness')
d_brain = subset(d_brain,select=extract_columns)#ICV
head(d_brain)

d2 = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/neuroCHARGE_WMH_ctxTH_age/WMH_ctxTH/Example_SYS/input_WMH_ICV_data_2019-09-16.tsv')
names(d2)

d_brain = d_brain %>% left_join(subset(d2,select=c('uniqueID',"ICVmL")) %>% dplyr::rename(ICV=ICVmL))
d_brain = d_brain %>% dplyr::select(-derived_lefthemvolume,-derived_righthemvolume)
dir_to_write = '/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts'
std_brain_colnames = fread(file.path(dir_to_write,'BrainPhenoNames.txt'),header=F)$V1
if(all(std_brain_colnames %in% names(d_brain))){#TRUE
	write_tsv(d_brain,file.path(dir_to_write,"SPS_BrainPhenos.tsv"))
}else{
	stop('check the column names of the BrainPheno data file.\n')
}

d_diabetes=fread(file.path(sys_data_dir,'jean_shin_dataanalysis_completespssubjects_final_20231108_170937.txt'))
names(d_diabetes) = str_remove(names(d_diabetes),"SPS-SPSm07_rawdata.")
with(d_diabetes,table(Diabetes,DiabetesType))
d_diabetes$DiabetesType = as.numeric(d_diabetes$DiabetesType)
d_diabetes = d_diabetes %>% 
	mutate(T2D = ifelse(Diabetes == 1 & DiabetesType == 3, TRUE,NA)) %>%
	mutate(T2D = ifelse(Diabetes == 2 & DiabetesType == -888, FALSE, T2D))
with(d_diabetes,table(DiabetesType,T2D,useNA='a'))
with(d_diabetes,table(Diabetes,T2D,useNA='a'))

d3 = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/neuroCHARGE_WMH_ctxTH_age/WMH_ctxTH/Example_SYS/input_SYS_Cov_data_wi_DM.tsv')
head(d3)

## bmi
d_bmi = fread(file.path(sys_data_dir,'jean_shin_dataanalysis_completespssubjects_20220916_151019.txt'))
strs_to_remove = c('SPS-N06f_body_anthropometry.','SPS-MRI-AbnominalFat_20150616.','SPS-N06f_body_bioimpedance.')
for(str_to_remove in strs_to_remove) {
  names(d_bmi) = str_remove(names(d_bmi),str_to_remove)
}
head(d_bmi)

d_bmi = d_bmi %>% dplyr::select(uniqueID,bmi) %>%
  mutate(bmi = as.numeric(bmi)) %>%
  mutate(logbmi = log(bmi))

tmp = remove.outliers_grubbs(d_bmi,"logbmi")
tmp$counts.NA
head(tmp$cleaned_data)
tmp$cleaned_data %>% filter(is.na(logbmi))
d_bmi = tmp$cleaned_data %>% 
  mutate(bmi = ifelse(is.na(logbmi),NA, bmi))
d_bmi %>% filter(is.na(logbmi))

## age
d_age = fread(file.path(sys_data_dir,'jean_shin_dataanalysis_completespssubjects_20220916_150904.txt'))
psych::describe(d_age)
names(d_age) = str_remove(names(d_age),'SPS-FieldModule.')
d_age %>% ggplot(aes(x=AgeParentAtHospitalVisitInMonths,y=AgeAtMRI_InMonths)) +
  geom_point() + 
  geom_abline(slope=1,intercept = 0)
head(d_age %>% mutate(time.reverse = AgeParentAtHospitalVisitInMonths-AgeAtMRI_InMonths) %>% arrange(time.reverse),20)
d_age = d_age %>% mutate(age_MRI = AgeAtMRI_InMonths/12, age_BMI = AgeParentAtHospitalVisitInMonths/12)
head(d_age)

## education
d_edu = fread(file.path(sys_data_dir,'jean_shin_dataanalysis_completesyssubjects_rev1_20231109_145908.txt'))
strs_to_remove = c('ADO1-A07f.','SPS-SYSParentsMerge.')
for (str_to_remove in strs_to_remove){
  names(d_edu) = str_remove(names(d_edu),str_to_remove)
}
# d_edu = d_edu %>% filter(Education>0)
# 1=primary not completed 
# 2=primary completed 
# 3=high school not completed 
# 4=high school completed 
# 5=college not completed 
# 6=college completed 
# 7=university not completed 
# 8=bachelor completed 
# 9=master or doctorate 
# 10=unknown (C1a)
head(d_edu)
table(d_edu$Education,useNA = 'a')
d_edu  = d_edu %>% 
  mutate(education_binary = ifelse(Education<0 | Education==10,NA,"High")) %>% 
  #high: higher than high-school completion
  mutate(education_binary = ifelse(Education >=1 & Education <=4,"Low",education_binary)) %>%
  mutate(education_binary = ifelse(Education >4 & Education <=9,"High",education_binary)) 
with(d_edu,table(Education,education_binary,useNA = 'a'))
colSums(with(d_edu,table(Education,education_binary)))
prop.table(colSums(with(d_edu,table(Education,education_binary))))
#       High          Low 
#419 (45.3%)  502 (54.7%) 

# merging datasets ------------------------------------------------------------
cov_to_anal = merge(subset(d3,select=c(uniqueID,age_years,Sex,current_smoking_Wave2SmokingStatus2,hypertension)),
                    subset(d2,select=c(uniqueID,ICVmL,Brain)),
                    by="uniqueID") %>%
              dplyr::select(-ICVmL,-Brain)%>%
              dplyr::rename(current_smoking = current_smoking_Wave2SmokingStatus2) %>% 
              left_join(d_diabetes %>% dplyr::select(uniqueID,T2D)) %>%
              filter(uniqueID %in% d_brain$uniqueID)              
cov_to_anal = cov_to_anal %>% left_join(d_bmi %>% dplyr::select(uniqueID,bmi))
cov_to_anal = cov_to_anal %>% left_join(d_age %>% dplyr::select(uniqueID,age_MRI,age_BMI))#age_years is the age_MRI
cov_to_anal = cov_to_anal %>% left_join(d_edu %>% dplyr::select(uniqueID,education_binary))#age_years is the age_MRI

head(cov_to_anal)#611 with MRI-age
psych::describe(cov_to_anal)
print(
  data.frame(var=rownames(psych::describe(cov_to_anal)),
             n=psych::describe(cov_to_anal)[['n']])
)

# write out the BMI_covariate file --------------------------------------------
write_tsv(cov_to_anal,file.path(dir_to_write,"SPS_BMI_covariates.tsv"))

# APOE4 and genotype PCs ------------------------------------------------------
APOE4 = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/Grant_NIH_Metabolon_ApoE_ZP/NIH_metabolon_ApoE/SYS_all_ApoEgeno.tsv')
APOE4 = APOE4 %>%
	mutate(E4_status = ifelse(str_detect(ApoE,"e4"),"E4-carrier","E4-noncarrier" )) %>%
	mutate(E4_status = ifelse(ApoE=="e2e4",NA,E4_status))

with(APOE4,table(ApoE,E4_status,useNA='a'))

APOE4 = APOE4 %>% right_join(d_brain %>% dplyr::select(uniqueID))
dim(APOE4)

## genotype PCs
genotypePCs = fread(file.path(sys_data_dir,'jean_shin_dataanalysis_completespssubjects_final_20231109_111838.txt'))
head(genotypePCs)
names(genotypePCs) = str_remove(names(genotypePCs),"SPS-PrincipalComponents.")
APOE4_wi_genoPCs = APOE4 %>% left_join(genotypePCs %>% dplyr::select(-fam_id,-sub_id))
head(APOE4_wi_genoPCs)

## write-out the files
write_tsv(APOE4_wi_genoPCs,file.path(dir_to_write,"SPS_APOE4_wi_genoPCs.txt"))
