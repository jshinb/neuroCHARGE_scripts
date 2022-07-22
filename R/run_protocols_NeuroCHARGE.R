# WMH-ctxTH NeuroCHARGE project ----------------------------------------------#
#
# Date Created: July 22, 2022
# Author: Jean Shin (jean.shin@sickkids.ca)
#
#
# ******************** SPECIFY INPUT FILE HERE ******************** 
#
#
input_specification_file = 
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

pacman::p_load(tidyverse, data.table, psych, readxl,
               patchwork, GGally,corrplot,hrbrthemes,
               tableone,FactoMineR,factoextra,
               ppcor, lsmeans, multcomp, ModelMetrics,
               caret, gridExtra, Hmisc, pastecs, testit)


#source the input-specification file ------------------------------------------
# need to think about what to do
op <- options(nwarnings = 10000)
setwd("./") #current directory
source(input_specification_file)

outdir=paste(cohort_name,ancestry,"outputs",sep="_")
dir.create(outdir)
file.copy(input_specification_file,outdir)

# source functions ------------------------------------------------------------
cat("Prep: sourcing functions\n")
source("0_functions_get_adj_ctxTH_WMH.R")

# step1: ----------------------------------------------------------------
source("1_DataWrangle.R")

# step2: ----------------------------------------------------------------------
source("2_Descriptive_statistics.R")

# step3: ----------------------------------------------------------------------
source("3a_roiAssociation_statistics_baseModel.R")
source("3b_roiAssociation_statistics_fullModel.R")
source("3c_plot_roiAssociation_results_base_vs_fullModels.R")

# step4: ----------------------------------------------------------------------
source("4_PC1_WMH_insulaTH.R")

# step5: ----------------------------------------------------------------------
source("5_roiAssociation_statistics_baseModel_AgeDeciles.R")

cat("\n# -------------------------------------------------------------------------------------- #\n",
    file=file.path(outdir,input_specification_file), append=T)
cat("Warnings:\n",file=file.path(outdir,input_specification_file), append=T)
capture.output(summary(warnings()),
               file=file.path(outdir,input_specification_file), append=T)
