setwd('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts')

# BMI-ctxVolume NeuroCHARGE project ------------------------------------------### Date Created: Nov 09, 2023# Author: Jean Shin (jean.shin@sickkids.ca)# SPS data# To-Do: add BMI-stratified [need to see sample sizes] ** UKBB# # ******************** SPECIFY INPUT FILE HERE ******************************** ##input_specification_file = 'cohort_specific_inputs_UKBB.txt'### *****************************************************************************#-----------------------------------------------------------------------------## install/load libraries ------------------------------------------------------cat("Prep: installing/loading libraries\n")# install pacman package for 'installing' and 'loading' packagesif(!is.element("pacman",installed.packages()[,1])){  install.packages("pacman")}pacman::p_load(tidyverse, tidyselect, data.table, readxl, psych, outliers,               patchwork, ggcorrplot, hrbrthemes, gridExtra, ggseg,               tableone, FactoMineR, factoextra,               ppcor, lsmeans, multcomp, ModelMetrics,               # GGally,caret,                Hmisc, pastecs, testit,               futile.logger,tryCatchLog)# for error/warning messagesif(length(find.package("ggsegDefaultExtra", quiet=TRUE, verbose=TRUE))==0){  remotes::install_github("LCBC-UiO/ggsegDefaultExtra")}library(ggsegDefaultExtra)#source the input-specification file ------------------------------------------# need to think about what to doop <- options(nwarnings = 10000, warn=0)setwd("./") #current directorysource(input_specification_file)outdir=paste(cohort_name,ancestry,"outputs",sep="_")# # remove the output directory if a previous folder exists# if(file.exists(outdir)){#   unlink(outdir,recursive = T)# }dir.create(outdir)file.copy(input_specification_file,outdir,overwrite = T)# define functions ------------------------------------------------------------cat("Prep: sourcing functions\n")source("0_functions_get_adj_ctxVol_BMI.R")df_filter = read_xlsx('df_filter.xlsx')# 1. data wrangling -----------------------------------------------------------logfile=file.path(outdir,"1_DataWrangle.log"); op = options(warn = 1)options(keep.source = TRUE)        # source code file name and line number trackingflog.appender(appender.file(logfile))  # to log into a file instead of consoleflog.threshold(INFO)    # TRACE, DEBUG, INFO, WARN, ERROR, FATALtryCatchLog(source("1_DataWrangle.R"))options(op)

# tableOne --------------------------------------------------------------------
# https://ehsanx.github.io/intro2R/data-summary-with-tableone.html
filter_i0 = 12; filter_i1 = 17
analdat = c()for(filter_i in c(filter_i0:filter_i1)){    group_name = df_filter$group_name[filter_i]    filter_expression_i = df_filter$expression[filter_i]	analdati = get_subset(data = "d",filter_expression=filter_expression_i) %>% 
	dplyr::select(IID,BMI,age,sex) %>%
	mutate(group = group_name)	
	analdat = rbind(analdat,analdati)
}
#dput(names(table(analdat$group)))

bmi.levels= c("BMI <20kg/m2", 
"20<= BMI <25kg/m2", 
"25<= BMI <30kg/m2", 
"30<= BMI <35kg/m2", 
"35<= BMI <40kg/m2", 
"BMI >=40kg/m2")

no_strata = CreateTableOne(data = analdat %>% mutate(group=factor(group,levels=bmi.levels)),
                         vars = c("age", "sex", "group"), 
                         factorVars = c("sex","group")
                         )

print(no_strata,
      showAllLevels = TRUE,
      nonnormal = "Age"
     )
tab_csv <- print(no_strata,
                 nonnormal = "age",
                 printToggle = FALSE)
write.csv(tab_csv, file = file.path(outdir,"tableOne_summary_age_sex_BMI_group.csv"))


strata <- CreateTableOne(data = analdat,
                         vars = c("age", "sex"), ## Note that BMI-group is not included because we already have strata = Gender
                         factorVars = c("sex"), ## Again, BMI-group is not included because it is in the strata argument
                         strata = "group" ## BMI-group
                         )
print(strata, nonnormal="age",cramVars = 'group')

tab_csv <- print(strata,
                 nonnormal = "age",
                 printToggle = FALSE)
write.csv(tab_csv, file = file.path(outdir,"tableOne_summary_age_sex_by_BMI_group.csv"))