#----------------------------------------------------
rm(list=ls());gc(T)

library(tidyverse)
library(data.table)

wd = "/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts"
setwd(wd)

f="/Users/jshin/Library/CloudStorage/OneDrive-SickKids/SYS_Database/jean_shin_dataanalysis_completesyssubjects_rev1_20221003_163611.txt"
d = fread(f)

strs_to_remove = c("ADO1-MRI-FreeSurferMay2012.",
                   "ADO1-MRI-AbdoFat_2012Oct23.",'
                   ADO1-N06f_body_anthropometry.',
                   'ADO1-N06f_body_bioimpedance.')
for(str_to_remove in strs_to_remove) {
  names(d) = str_remove(names(d),str_to_remove)
}
rois = c("bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", 
         "cuneus", "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal", 
         "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal", 
         "lingual", "medialorbitofrontal", "middletemporal", "parahippocampal", 
         "paracentral", "parsopercularis", "parsorbitalis", "parstriangularis", 
         "pericalcarine", "postcentral", "posteriorcingulate", "precentral", 
         "precuneus", "rostralanteriorcingulate", "rostralmiddlefrontal", 
         "superiorfrontal", "superiorparietal", "superiortemporal", "supramarginal", 
         "frontalpole", "temporalpole", "transversetemporal", "insula")  
"LhbanksstsVolume"
vol_names = paste("Lh",rois,"Volume",sep='')
vol_names = c(vol_names,paste("Rh",rois,"Volume",sep=''))
length(intersect(names(d),vol_names))#68

area_names = paste("Lh",rois,"Area",sep='')
area_names = c(area_names,paste("Rh",rois,"Area",sep=''))
length(intersect(names(d),area_names))#68

th_names = paste("Lh",rois,"Thickness",sep='')
th_names = c(th_names,paste("Rh",rois,"Thickness",sep=''))
th_names = c(th_names,paste("Thickness",c("Lh","Rh"),"emisph",sep=''))
length(intersect(names(d),th_names))#70

d_to_write = subset(d,select=c('uniqueID',"fam_id","QC",vol_names,area_names,th_names))
names(d_to_write)
d_to_write = data.frame(d_to_write)
replace_neg_or_failedQC_wi_na = function(x,QC){
  ind = !is.na(x) & x<0
  ind = ind | (!is.na(QC) & QC==1)
  x[ind] = NA
  x
}

d_to_write[,c(vol_names,area_names,th_names)] <- apply(d_to_write[,c(vol_names,area_names,th_names)],2,as.numeric)
d_to_write[,c(vol_names,area_names,th_names)] <- apply(d_to_write[,c(vol_names,area_names,th_names)],2,replace_neg_or_failedQC_wi_na,QC=d_to_write$QC)
d_to_write %>% filter(QC==1)
psych::describe(d_to_write) %>% dplyr::select(n,min,max) %>% arrange(n)#988~989

rm_sample_ind = apply(apply(d_to_write[,c(vol_names,area_names,th_names)],2,is.na),1,sum) == ncol(d_to_write[,c(vol_names,area_names,th_names)])
table(rm_sample_ind )#40 individuals with no brain information

d_to_write = d_to_write %>% filter(!rm_sample_ind)#989

setdiff (names(d),names(d_to_write))

#BMI_Z_Score, ICV
# ICV: /Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/data_to_send/adolescent_brain_adiposity_data_2022-10-07.tsv
d_other = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/data_to_send/adolescent_brain_adiposity_data_2022-10-07.tsv')
d_BMI_cov = d_ICV %>% dplyr::select(uniqueID,age,sex,BMI_Z_Score)

