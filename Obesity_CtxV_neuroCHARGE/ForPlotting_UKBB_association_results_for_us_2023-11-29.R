rm(list=ls());gc(T)
setwd("/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts")
library(tidyverse)
library(data.table)


fnames='roi_assoc_res_two_methods_base_old_heavy.Rdata
roi_assoc_res_two_methods_base_old_light.Rdata
roi_assoc_res_two_methods_base_young_heavy.Rdata
roi_assoc_res_two_methods_base_young_light.Rdata
roi_assoc_res_two_methods_full_old_heavy.Rdata
roi_assoc_res_two_methods_full_old_light.Rdata
roi_assoc_res_two_methods_full_young_heavy.Rdata
roi_assoc_res_two_methods_full_young_light.Rdata'

fnames=unlist(str_split(fnames,'\n'))
mod_subset_names = str_remove(fnames,"roi_assoc_res_two_methods_")
mod_subset_names = str_remove(mod_subset_names,".Rdata")

res1_all = c()
for(fi in 1:length(fnames)){
  load(file.path("UKBB2_EUR_outputs",fnames[fi]))
  head(res1)
  res1 = res1 %>% mutate(mod=unlist(str_split(mod_subset_names[fi],"_"))[1])
  res1 = res1 %>% mutate(group = paste0(unlist(str_split(mod_subset_names[fi],"_"))[-1],collapse ="_"))
  res1_all = rbind(res1_all,res1)
  rm(res1)
}

res1_all = res1_all %>% filter(!str_detect(E4_status,'interaction'))
head(res1_all)

res1_all = res1_all %>% mutate(U95CI = Estimate + 1.96*SE,L95CI = Estimate - 1.96*SE)
range(res1_all$L95CI,res1_all$U95CI)#
