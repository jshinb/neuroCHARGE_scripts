# BMI-ctxVolume NeuroCHARGE project ------------------------------------------#
#
# Date Created: Nov 09, 2023
# Author: Jean Shin (jean.shin@sickkids.ca)
# SPS data
# To-Do: add BMI-stratified [need to see sample sizes] ** UKBB
# 
# ******************** SPECIFY INPUT FILE HERE ******************************** 
#
#
input_specification_file = 'cohort_specific_inputs_ukb677610.txt'
#
#
# *****************************************************************************
#-----------------------------------------------------------------------------#

# install/load libraries ------------------------------------------------------
cat("Prep: installing/loading libraries\n")

# install pacman package for 'installing' and 'loading' packages
if(!is.element("pacman",installed.packages()[,1])){
  install.packages("pacman")
}

pacman::p_load(tidyverse, tidyselect, data.table, readxl, psych, outliers,
               patchwork, ggcorrplot, hrbrthemes, gridExtra, ggseg,
               tableone, FactoMineR, factoextra,
               ppcor, lsmeans, multcomp, ModelMetrics,
               # GGally,caret, 
               Hmisc, pastecs, testit,
               futile.logger,tryCatchLog)# for error/warning messages


if(length(find.package("ggsegDefaultExtra", quiet=TRUE, verbose=TRUE))==0){
  remotes::install_github("LCBC-UiO/ggsegDefaultExtra")
}

library(ggsegDefaultExtra)

#source the input-specification file ------------------------------------------
# need to think about what to do
op <- options(nwarnings = 10000, warn=0)
setwd("./") #current directory
source(input_specification_file)

outdir=paste(cohort_name,ancestry,"outputs",sep="_")
# # remove the output directory if a previous folder exists
# if(file.exists(outdir)){
#   unlink(outdir,recursive = T)
# }
dir.create(outdir)
file.copy(input_specification_file,outdir,overwrite = T)

# define functions ------------------------------------------------------------
cat("Prep: sourcing functions\n")
source("0_functions_get_adj_ctxVol_BMI.R")
df_filter = read_xlsx('df_filter.xlsx')

# 1. data wrangling -----------------------------------------------------------
logfile=file.path(outdir,"1_DataWrangle.log"); 
op = options(warn = 1)
options(keep.source = TRUE)        # source code file name and line number tracking

flog.appender(appender.file(logfile))  # to log into a file instead of console
flog.threshold(INFO)    # TRACE, DEBUG, INFO, WARN, ERROR, FATAL

tryCatchLog(source("1_DataWrangle.R"))
options(op)

# 2. descriptive stats --------------------------------------------------------
d = d %>% dplyr::select(-matches("temporalpole"))
rois = c("bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", 
         "cuneus", "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal", 
         "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal", 
         "lingual", "medialorbitofrontal", "middletemporal", "parahippocampal", 
         "paracentral", "parsopercularis", "parsorbitalis", "parstriangularis", 
         "pericalcarine", "postcentral", "posteriorcingulate", "precentral", 
         "precuneus", "rostralanteriorcingulate", "rostralmiddlefrontal", 
         "superiorfrontal", "superiorparietal", "superiortemporal", "supramarginal", 
         "frontalpole", "transversetemporal", "insula")  

source("2_Descriptive_statistics.R")
#ukb677610 - does not have FreeSurfer information on temporalpole
#define data groupings for subset association analyses
df_filter = read_xlsx('df_filter.xlsx')
getN = function(data_name, filter_expression){
  execute_string = paste("analdat <- ",data_name,"%>% filter(",filter_expression,")",sep='')
  eval(parse(text=execute_string))
  analdat
  ret = nrow(analdat)
  ret
}

NE4_carrier <- nrow(d %>% filter(E4_status=="E4-carrier"))
NE4_noncarrier <- nrow(d %>% filter(E4_status=="E4-noncarrier"))
filter_is = 2:nrow(df_filter)
for(filter_i in filter_is){
  filter_expression_i = df_filter$expression[filter_i]
  analdat = get_subset("d",filter_expression = filter_expression_i)
  NE4_carrier = c(NE4_carrier,nrow(analdat %>% filter(E4_status=="E4-carrier")))
  NE4_noncarrier = c(NE4_noncarrier,nrow(analdat %>% filter(E4_status=="E4-noncarrier")))
  rm(analdat)
}

df_filter = df_filter %>% mutate(N_E4_noncarrier=NE4_noncarrier,N_E4_carrier=NE4_carrier,index=1:nrow(df_filter))
cat(df_filter$group_name,sep='\n')

# 3. regional associations ----------------------------------------------------
source("3_roiAssociation_statistics_group_model.R")
source("3_plot_roiAssociation_results_base_vs_fullModels.R")

###* REMOVE THIS LINE BEFORE DISTRIBUTION:
###* need to do this after looking at the results from other cohorts
# 4. PCA of BMI and insular thickness adjusted for basic covariates -----------
source("4_PC1_BMI_insulaTH.R")
source("4b_PC1_BMI_MCT.R")
