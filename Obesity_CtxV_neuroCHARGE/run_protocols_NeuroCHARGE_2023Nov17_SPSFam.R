# BMI-ctxVolume NeuroCHARGE project ------------------------------------------#
#
# Date Created: Nov 09, 2023
# Author: Jean Shin (jean.shin@sickkids.ca)
# 
# ******************** SPECIFY INPUT FILE HERE ******************************** 
#
#
input_specification_file = 'cohort_specific_inputs_SPSFam.txt'
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
               futile.logger,tryCatchLog,
               lmerTest)# for error/warning messages


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
source("0_functions_get_adj_ctxVol_BMI_FAM.R")

# 1. data wrangling -----------------------------------------------------------
logfile=file.path(outdir,"1_DataWrangle.log"); 
op = options(warn = 1)
options(keep.source = TRUE)        # source code file name and line number tracking

flog.appender(appender.file(logfile))  # to log into a file instead of console
flog.threshold(INFO)    # TRACE, DEBUG, INFO, WARN, ERROR, FATAL

tryCatchLog(source("1_DataWrangle_FAM.R"))
options(op)

# 2. descriptive stats --------------------------------------------------------
source("2_Descriptive_statistics.R")

# 3. regional associations ----------------------------------------------------
source("3_roiAssociation_statistics_group_model_FAM.R")
source("3_plot_roiAssociation_results_base_vs_fullModels.R")

