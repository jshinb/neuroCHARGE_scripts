# Select individuals -------------------------------------------------------------------
# unrelated
# White British
# no history of stroke and dementia

# 
run.outsideR <- function(){

module load gcc/8.3.0 intel/2019u4 r/4.1.2;R

}

library(ukbtools)
library(stringr)
library(dplyr)
library(data.table)
library(readr)

setwd("/scratch/t/tpaus/jshinb/UKBB/datasets/ukb677610/")

# read in data ------------------------------------------------------------------------
d = ukb_df('ukb677610')
d_brain_variables = fread('/scratch/t/tpaus/jshinb/KEEP/UKBB_JS/desikan_white_variables.txt',header=F)$V1
d_brain_variables = d_brain_variables[str_detect(d_brain_variables,"_2_0")]
d_brain2 = d %>% dplyr::select(all_of(c("eid",d_brain_variables)))

# selectID 
stroke.fid = 'f42006_0_0'#date of stroke
AD.fid = 'f42020_0_0' #date of AD
all.cause.dementia.fid ='f42018_0_0'#date.date
white.british.fid = 'f22006_0_0'#409418 vs. 92793 non-british white
aneuploidy = 'f22019_0_0'
heterozygosity_or_missingness = 'f22027_0_0'
missing_call_rate = 'f22005_0_0'

selectID_strs = c(stroke.fid,AD.fid,all.cause.dementia.fid,white.british.fid,aneuploidy,heterozygosity_or_missingness,missing_call_rate)

selectID_cols_ind = 'eid'
for(i in selectID_strs){
selectID_cols_ind = c(selectID_cols_ind,names(d)[str_detect(names(d),i)])
}

d_selectID = d %>% dplyr::select(all_of(selectID_cols_ind)); n0 = nrow(d_selectID); 
print(n0)#502211

d_selectID = d_selectID %>% filter(genetic_ethnic_grouping_f22006_0_0==1); 
n1 = nrow(d_selectID); 
print(c(n1,n0-n1)); #92,793
n0=n1 #remove 92,793 non-british white

d_selectID = d_selectID %>% filter(is.na(date_of_all_cause_dementia_report_f42018_0_0))#remove 8313 with all-cause dimentia
n1 = nrow(d_selectID); 
print(c(n1,n0-n1)); # 8,313 removed
n0=n1 #401,105 left
head(d_selectID %>% filter(!is.na(date_of_alzheimers_disease_report_f42020_0_0)))#none

d_selectID = d_selectID %>% filter(is.na(date_of_stroke_f42006_0_0)) #remove 15156 with stroke
n1 = nrow(d_selectID); 
print(c(n1,n0-n1)); # 15,156 removed
n0=n1 #385,949 left

d_selectID = d_selectID %>% filter(is.na(sex_chromosome_aneuploidy_f22019_0_0)) #remove XX with aneuploidy
n1 = nrow(d_selectID); 
print(c(n1,n0-n1)); # 506 removed
n0=n1 #385,443 left

d_selectID = d_selectID %>% filter(is.na(outliers_for_heterozygosity_or_missing_rate_f22027_0_0)) #remove XX with heterozygosity or missingness
n1 = nrow(d_selectID); 
print(c(n1,n0-n1)); # 685 removed
n0=n1 #384,758 left
#summary(d_selectID %>% filter(is.na(outliers_for_heterozygosity_or_missing_rate_f22027_0_0)) %>% select(missingness_f22005_0_0))
 
d_genetic_sex = fread('/home/t/tpaus/jshinb/ukbb_gwas_data/ukbb41449_genetic_sex_2023-10-20.tsv')
d_report_sex = d %>% select(all_of(c("eid","sex_f31_0_0"))) %>% 
	mutate(report_sex = ifelse(sex_f31_0_0==0, 'Female', "Male"))#0 female
d_report_sex = d_report_sex %>% left_join(d_genetic_sex)
d_selectID = d_selectID %>% left_join(d_report_sex)#n=385,949
d_selectID = d_selectID %>% filter(report_sex==genetic_sex_f22001_0_0)#385,649 after removing 300 individuals
n1 = nrow(d_selectID); 
print(c(n1,n0-n1)); # 163 removed
n0=n1 #384,595 left

#read in APOE4 data
d_apoE = fread('/scratch/t/tpaus/jshinb/KEEP/UKBB_JS/hardcalled_genotypes_for_APOE_genotypes_2023-05-04.tsv') %>%
  dplyr::rename(eid=IID)#44,890

# read in the brain data
d_brain = fread('/scratch/t/tpaus/jshinb/KEEP/UKBB_JS/ukb677610_d_desikan_white_variables_2024-04-02.txt')
d_brain = d_brain %>% left_join(d_selectID %>% dplyr::select(eid)) %>% filter(eid %in% d_apoE$eid)
n1 = nrow(d_brain)
print(c(n1,n0-n1)); # 339,062 without brain and/or APOE data removed
n0=n1 #45,533 left

# remove related individuals -----------------------------------------------------------
my_relatedness_data = fread('/project/t/tpaus/tpaus/UKBB/datasets/gene_data/ukb43688_rel_s488264.dat')
rmID = ukb_gen_samples_to_remove(my_relatedness_data, d_brain$eid, cutoff = 0.0884)#643 individuals
d_brain_unrelated = d_brain %>% filter(!eid %in% rmID)
n1 = nrow(d_brain_unrelated)
print(c(n1,n0-n1)); # 643 related individuals removed
n0=n1 #44,890 left
d_brain_unrelated2 = d_brain_unrelated %>% dplyr::select(eid) %>% left_join(d_brain2)

my_relatedness_data %>% filter(ID1 %in% d_brain$eid,ID2 %in% d_brain$eid)
exp.ID = unique(unlist(c(my_relatedness_data %>% filter(ID1=='1880450'|ID2=="1880450") %>% select(ID1,ID2))))
d_brain %>% filter(eid %in% exp.ID) %>% select(eid)
d_brain_unrelated %>% filter(eid %in% exp.ID) %>% select(eid) #5303405
# add icv
d_icv = d %>% dplyr::select(matches("eid|f26521"))
d_brain_unrelated = d_brain_unrelated %>% left_join(d_icv)
d_brain_unrelated2 = d_brain_unrelated2 %>% left_join(d_icv)
d_brain = d_brain2; rm(d_brain2)
d_brain_unrelated = d_brain_unrelated2; rm(d_brain_unrelated2)

d_column_names = str_split(names(d_brain_unrelated),"_f",simplify=T)
d_column_names = data.frame(d_column_names)
d_column_names %>% filter(X3!="")
d_column_names = d_column_names %>% mutate(X4=ifelse(X3!='',paste(X1,"_f",X2,sep=''),X1))
d_column_names = d_column_names %>% mutate(X5 = str_replace(X4,"left_hemisphere","lh"))
d_column_names = d_column_names %>% mutate(X5 = str_replace(X5,"right_hemisphere","rh"))
d_column_names = d_column_names %>% mutate(X5 = str_replace(X5,"volume_of","volume"))
d_column_names = d_column_names %>% mutate(X5 = str_replace(X5,"area_of","area"))
d_column_names = d_column_names %>% mutate(X5 = str_replace(X5,"mean_thickness_of","thickness"))
d_column_names$X5 = str_replace(d_column_names$X5,"whole_brain","wholebrain")
d_column_names$X5 = str_replace(d_column_names$X5,"thickness_globalmeanmean","globalmeanmean")
d_column_names = cbind(d_column_names,str_split(d_column_names$X5,"_",simplify=T))
d_column_names = data.frame(d_column_names)
d_column_names = d_column_names %>% mutate(col_names = ifelse(X3.1!="",paste(X3.1,X2.1,X1.1,sep="_"),X1))

ind = str_detect(names(d_brain_unrelated),"_insula_")
ind = ind &str_detect(names(d_brain_unrelated),"volume")
ind = ind &str_detect(names(d_brain_unrelated),"_2_0")
names(d_brain_unrelated)[ind]
names(d_brain_unrelated) = d_column_names$col_names
names(d_brain_unrelated)[ind]
readr::write_tsv(d_brain_unrelated,"/scratch/t/tpaus/jshinb/KEEP/UKBB_JS/ukb677610_d_brain_unrelated.txt")
