# BMI-ctxVolume NeuroCHARGE project ------------------------------------------#
#
# Date Created: Oct 23, 2023
# Author: Jean Shin (jean.shin@sickkids.ca)
#
#
# ******************** SPECIFY INPUT FILE HERE ******************** 
#
#
input_specification_file = 'cohort_specific_inputs_UKBB.txt'
#
#
# *****************************************************************
#-----------------------------------------------------------------------------#

# install/load libraries ------------------------------------------------------
cat("Prep: installing/loading libraries\n")

# install pacman package for 'installing' and 'loading' packages
if(!is.element("pacman",installed.packages()[,1])){
  install.packages("pacman")
}

pacman::p_load(tidyverse, tidyselect, data.table, readxl, psych, outlier,
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

# remove the output directory if a previous folder exists
if(file.exists(outdir)){
  unlink(outdir,recursive = T)
}
dir.create(outdir)
file.copy(input_specification_file,outdir,overwrite = T)

# define functions ------------------------------------------------------------
cat("Prep: sourcing functions\n")
source("0_functions_get_adj_ctxVol_BMI.R")

# 1. data wrangling -----------------------------------------------------------
logfile=file.path(outdir,"1_DataWrangle.log"); 
op = options(warn = 1)
options(keep.source = TRUE)        # source code file name and line number tracking

flog.appender(appender.file(logfile))  # to log into a file instead of console
flog.threshold(INFO)    # TRACE, DEBUG, INFO, WARN, ERROR, FATAL

tryCatchLog(source("1_DataWrangle.R"))
options(op)

# 2. descriptive stats --------------------------------------------------------
source("2_Descriptive_statistics.R")

# 3. regional associations ----------------------------------------------------
source("3a_roiAssociation_statistics_baseModel.R")
source("3b_roiAssociation_statistics_fullModel.R")
source("3c_plot_roiAssociation_results_base_vs_fullModels.R")

source("3a_age_lt_65_roiAssociation_statistics_baseModel.R")
source("3b_age_lt_65_roiAssociation_statistics_fullModel.R")
source("3c_age_lt_65_plot_roiAssociation_results_base_vs_fullModels.R")

source("3a_age_gt_65_roiAssociation_statistics_baseModel.R")
source("3b_age_gt_65_roiAssociation_statistics_fullModel.R")
source("3c_age_gt_65_plot_roiAssociation_results_base_vs_fullModels.R")

###* REMOVE THIS LINE BEFORE DISTRIBUTION:
###* need to do this after looking at the results from other cohorts
# 4. PCA of BMI and insular thickness adjusted for basic covariates -----------
source("4_PC1_BMI_insulaTH.R")
source("4b_PC1_BMI_MCT.R")
