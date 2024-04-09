#' ---
#' title: "sBMI-ctxVolume NeuroCHARGE project"
#' author: "Jean Shin"
#' date created: "Nov 09, 2023"
#' date: "Jan 25, 2024"
#' ---

# To-Do: add BMI-stratified [need to see sample sizes] ** UKBB
# 
# ******************** SPECIFY INPUT FILE HERE ******************************** 
#
#
input_specification_file = 'cohort_specific_inputs_UKBB.txt'
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

# fitting ---------------------------------------------------------------------
#https://cran.r-project.org/web/packages/tidymv/vignettes/plot-smooths.html
theme_set(theme_bw())
library(mgcv)
library(tidygam)
library(tidymv)

fit1 = gam(insula_thickness ~ s(age,k=4, by=E4_status) + s(age, k=4,by=sex)+ 
             s(BMI, k=4,by=E4_status) + s(BMI, k=4,by=sex)+education + mri_site, data=d)
par(mfrow=c(3,3));plot(fit1)

# plotting smooth curves ------------------------------------------------------
plot = FALSE
d = d %>% mutate(sex = factor(sex,levels=c("M","F")),
                 E4_status = factor(E4_status,levels=c("E4-noncarrier","E4-carrier")))
d = d %>% mutate(sex.E4 = interaction(sex,E4_status))
pheno = 'thickness'
p_list = vector(mode='list',length=length(rois));i=0
for(region in rois){
  i=i+1
  region_pheno = paste(region,pheno,sep="_")
  #s(age, k=4,by=sex)+s(BMI, k=4,by=sex) + 
  mod = paste(region_pheno, "s(age.jitter) + s(BMI) + s(age.jitter,k=4, by=E4_status) +  
             s(BMI, k=4,by=E4_status) + sex + E4_status + education + mri_site",sep="~")
  fit1 = gam(as.formula(mod), data=d %>% mutate(age.jitter = jitter(age)))
  # mod.bam = "sex.E4 +s(BMI,k=4,by=sex.E4)"
  mod.bam = "sex.E4 +s(BMI,k=4,by=sex.E4) +s(age,k=4,by=sex.E4)"
  mod.bam = as.formula(paste(region_pheno,mod.bam,sep="~"))
  fit1 = bam(as.formula(mod.bam), data=d)
  
  p_age = plot_smooths(
    model = fit1,
    series = age,
    comparison = E4_status,
    facet_terms = sex,
    split = list(sex.E4=c("sex","E4_status"))
  )+ scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
    scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
    theme(legend.position = "top") + ggtitle(region_pheno)
  p_bmi = plot_smooths(
    model = fit1,
    series = BMI,
    comparison = E4_status,
    facet_terms = sex,
    split = list(sex.E4=c("sex","E4_status"))
  ) + scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
    scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
    theme(legend.position = "top") #+ ggtitle(region_pheno)
  p_list[[i]] = p_bmi
  
  #printing
  if(plot){
    p = p_age+p_bmi+plot_layout(nrow=2, guides="collect")
    pfile=paste("gam_plots_",region_pheno,"_",Sys.Date(),".png",sep='')
    png(pfile,width=11,height=8.5,units='in',dpi=300)
    print(p)
    dev.off()
  }
}
i=4;p_list[[3*i-2]]/p_list[[3*i-1]]/p_list[[3*i]]
#render to produce html (without RMD file) ------------------------------------
#rmarkdown::render("explore_obesity_Cth_association_2024-01-25.R")