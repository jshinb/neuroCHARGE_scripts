# WMH-ctxTH NeuroCHARGE project -----------------------------------------------
# Date: April 27, 2022
# Author: Jean Shin (jean.shin@sickkids.ca)
# Note: some scripts are adapted from ENIGMA Brain age GWAS project
# 2. Association tests between WMH and cortical thickness across 34 regions
#
#
# ****SPECIFY INPUT FILE HERE:
#
input_specification_file=
#
#-----------------------------------------------------------------------------#

#source the input-specification file ------------------------------------------
# need to think about what to do
op <- options(nwarnings = 10000)
setwd("./") #current directory
source(input_specification_file)

outdir=paste(cohort_name,ancestry,"outputs",sep="_")
dir.create(outdir)
file.copy(input_specification_file,outdir)

#create log file --------------------------------------------------------------
messages=file(file.path(outdir,"all_messages.log"), open="wt")
sink(messages, type="message")
sink(messages, type="output")


# install/load libraries ------------------------------------------------------
cat("Prep: installing/loading libraries\n")

# install pacman package for 'installing' and 'loading' packages
if(!is.element("pacman",installed.packages()[,1])){
  install.packages("pacman")
}

pacman::p_load(tidyverse, data.table, psych, 
               patchwork, GGally,corrplot,hrbrthemes,
               tableone,FactoMineR,factoextra,
               ppcor, lsmeans, multcomp, ModelMetrics,
               caret, gridExtra, Hmisc, pastecs, testit)


# source functions ------------------------------------------------------------
cat("Prep: sourcing functions\n")
source("0_functions_get_adj_ctxTH_WMH.R")

# step1: ----------------------------------------------------------------
source("1_DataWrangle.R")

# step2: ----------------------------------------------------------------------
source("2_Descriptive_statistics.R")

# step3: ----------------------------------------------------------------------
source("3_roiAssociation_statistics.R")

# step4: ----------------------------------------------------------------------
source("4_PC1_WMH_insulaTH.R")

cat("\n# -------------------------------------------------------------------------------------- #\n",
    file=file.path(outdir,input_specification_file), append=T)
cat("Warnings:\n",file=file.path(outdir,input_specification_file), append=T)
capture.output(summary(warnings()),
               file=file.path(outdir,input_specification_file), append=T)
sink()
closeAllConnections()
print(readLines(file.path(outdir,"all_messages.log")))
