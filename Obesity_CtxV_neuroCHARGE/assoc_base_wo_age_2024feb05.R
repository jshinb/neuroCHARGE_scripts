## ----setup, include=FALSE---------------------------------------------------------------------------
# had some issues: https://support.posit.co/hc/en-us/articles/200534577-Resetting-RStudio-Desktop-s-State
#date created: Nov 09, 2023
# markdown examples
#http://www.gavindouglas.ca/pages/R_notebook_example.html#Highlighted_examples
# output: 
#   html_notebook:
#     code_folding: show
#     theme: cerulean
#     toc: true
#     toc_float:
#       collapsed: true

knitr::opts_chunk$set(echo = TRUE, root.dir="./")


## ----install_load_packages, include=TRUE------------------------------------------------------------
input_specification_file = 'cohort_specific_inputs_UKBB.txt'
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


## ----read_inputs------------------------------------------------------------------------------------
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


## ----define_functions-------------------------------------------------------------------------------

cat("Prep: sourcing functions\n")
source("0_functions_get_adj_ctxVol_BMI.R")



## ----data_wrangle-----------------------------------------------------------------------------------

logfile=file.path(outdir,"1_DataWrangle.log"); 
op = options(warn = 1)
options(keep.source = TRUE)        # source code file name and line number tracking

flog.appender(appender.file(logfile))  # to log into a file instead of console
flog.threshold(INFO)    # TRACE, DEBUG, INFO, WARN, ERROR, FATAL

tryCatchLog(source("1_DataWrangle.R"))
options(op)



## ----descriptives-----------------------------------------------------------------------------------

source("2_Descriptive_statistics.R")



## ----subset_categories------------------------------------------------------------------------------

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

outdir2 = "adj_base_wi_group_adj_cov_before_splitting_no_age_2024-02-05"
dir.create(file.path(outdir,outdir2))
adj_before_splitting = TRUE
do.smooth.plot = FALSE

## ----assoc_tests_adj_base,eval=TRUE-----------------------------------------------------------------

get_adj_ctx <- 
  function(roii_pheno,d,covariate_names){
    d = d %>% mutate(age.c = scale(age)[,1],
                     age.c2 = (scale(age)/2)^2)# for stability (Gelman)
    analdati = subset(d,select=c(covariate_names,roii_pheno)) 
    
    ## remove outliers in the original variable
    y_outliers_removed0 = remove.outliers_grubbs(analdati,varnames=roii_pheno)
    print(y_outliers_removed0$counts.NA)
    analdati[[roii_pheno]] = y_outliers_removed0$cleaned_data[[roii_pheno]]
    analdati[['y']] = inormal(analdati[[roii_pheno]])# already scaled
    
    covariate_names = covariate_names[which(!psych::describe(subset(analdati,select=covariate_names))$sd==0)]
    mod0.roii = paste0(c(setdiff(covariate_names,c('age','sex'))),
                       collapse = " + ")
    if(length(unique(analdati$sex))==2){
      ## adjusting for age, sex and other covariates:
      mod0.roii = paste("y ~ sex*(",mod0.roii,")",sep="")
      
    }else{#length(unique(analdati$sex))=1
      mod0.roii = paste("y ~ ",mod0.roii,sep="")
    }
    
    fit0.roii = lm(as.formula(mod0.roii),data=analdati,na.action=na.exclude)
    adj_yname=paste(roii_pheno,"adj",sep=".")
    analdati[[adj_yname]] = resid(fit0.roii)
    
    ## remove outliers in the adjusted variable 
    y_outliers_removed = remove.outliers_grubbs(analdati,varnames=adj_yname)
    analdati[[adj_yname]] = y_outliers_removed$cleaned_dat[[adj_yname]]
    ret = subset(analdati,select=adj_yname)
    
    cat("\nmodel: ",mod0.roii,"\n")
    cat("\nnumbers of NA\'s in the original and the cleaned ",roii_pheno," are:\n",sep='')
    print(y_outliers_removed$counts.NA)
    
    ret
  }

get_adj_BMI <- 
  function(d,covariate_names){
    analdati = subset(d,select=c('IID','BMI',covariate_names)) %>%
      # mutate(age.c = scale(age)[,1]/2) %>% mutate(age.c2=age.c^2) %>%
      mutate(logBMI = log(BMI))#for removing outliers
    
    # remove outliers based on bmi
    x_outliers_removed0 = remove.outliers_grubbs(analdati,varnames="logBMI")
    print(x_outliers_removed0$counts.NA)
    analdati[["logBMI"]] = x_outliers_removed0$cleaned_data[["logBMI"]]
    analdati[['x']] = inormal(analdati[["logBMI"]])
    covariate_names = covariate_names[which(!psych::describe(subset(analdati,select=covariate_names))$sd==0)]
    
    mod0 = paste0(c(setdiff(covariate_names,c('age','sex',"ICV"))),collapse = " + ")
    if(length(unique(analdati$sex))==2){
      mod0 = paste("x ~ sex*(",mod0,")",sep="")
    }else{
      mod0 = paste("x ~ ",mod0,sep="")
    }
    # age, (sex) and education
    adj_xname = paste("BMI","adj",sep=".")
    fit0 = lm(as.formula(mod0),data=analdati,na.action=na.exclude)
    analdati[[adj_xname]] = resid(fit0)
    
    ## remove outliers in the adjusted variable 
    x_outliers_removed = remove.outliers_grubbs(analdati,varnames=adj_xname)
    analdati[[adj_xname]] = x_outliers_removed$cleaned_dat[[adj_xname]]
    ret = subset(analdati,select=adj_xname)
    
    cat("\nmodel: ",mod0,"\n")
    cat("\nnumbers of NA\'s in the original and the cleaned BMI are:\n",sep='')
    print(x_outliers_removed$counts.NA)
    
    ret = subset(analdati,select=adj_xname)
    ret
  }
op <- options(warn=1)
rois = c("bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", 
         "cuneus", "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal", 
         "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal", 
         "lingual", "medialorbitofrontal", "middletemporal", "parahippocampal", 
         "paracentral", "parsopercularis", "parsorbitalis", "parstriangularis", 
         "pericalcarine", "postcentral", "posteriorcingulate", "precentral", 
         "precuneus", "rostralanteriorcingulate", "rostralmiddlefrontal", 
         "superiorfrontal", "superiorparietal", "superiortemporal", "supramarginal", 
         "frontalpole", "temporalpole", "transversetemporal", "insula")  

# starting --------------------------------------------------------------------# 
cov_additional= c('current_smoking','hypertension','type2diabetes')
roi_ctx_measurements = c(paste(rois,"volume",sep="_"),
                         paste(rois,"area",sep="_"),
                         paste(rois,"thickness",sep="_"))
#********************* NEED TO CHANGE BEFORE DISTRIBUTION *********************#
#
filter_is = df_filter$index
#
#******************************************************************************#
theme_set(theme_bw())
library(mgcv)
library(tidygam)
library(tidymv)

models = "base"
for(model in models){
  if(model=="base"){
    cov_BMI = setdiff(std_covariate_names,c("BMI","age_BMI","E4_status","ICV",cov_additional))
    cov_ctx = setdiff(std_covariate_names,c("BMI","age_BMI","E4_status",cov_additional))
  }else{#full
    cov_ctx = setdiff(std_covariate_names,c("BMI","age_BMI","E4_status"))
    cov_BMI = setdiff(std_covariate_names,c("BMI","age_BMI","E4_status","ICV"))
  }
  for(filter_i in filter_is){
    group_name = df_filter$group_name[filter_i]
    filter_group = df_filter$group_label[filter_i]
    filter_expression_i = df_filter$expression[filter_i]
    
    msg = paste("\n3. Start estimating associations between BMI vs. roi-ctx phenotype under ",
                model," model in the participants with ",group_name,".\n",sep='')
    cat(msg)
    logfile=file.path(outdir,paste("3_",filter_group,"_roiAssociation_statistics_",model,".log",sep='')); messages=file(logfile, open="wt")
    #sink(messages, type="message")
    #sink(messages, type="output")
    
    if(adj_before_splitting){
      adj_ctx = do.call(cbind,lapply(roi_ctx_measurements,get_adj_ctx,
                                     d=d,
                                     covariate_names = cov_ctx[!str_detect(cov_ctx,"age")]))
      adj_BMI = get_adj_BMI(d=d, covariate_names = cov_BMI[!str_detect(cov_BMI,"age")])
      
      analdat = data.frame( subset(d,select = c("IID",std_covariate_names)),adj_ctx,adj_BMI )
      
      analdat = get_subset("analdat",filter_expression = filter_expression_i)
    }else{
      analdat = get_subset("d",filter_expression = filter_expression_i)
      adj_ctx = do.call(cbind,lapply(roi_ctx_measurements,get_adj_ctx,
                                     d=analdat,
                                     covariate_names = cov_ctx[!str_detect(cov_ctx,"age")]))
      adj_BMI = get_adj_BMI(d=analdat, covariate_names = cov_BMI[!str_detect(cov_BMI,"age")])
      analdat = data.frame( subset(analdat,select = c("IID",std_covariate_names)),adj_ctx,adj_BMI )
    }
    
    if(min(table(analdat$E4_status)) >(30+length(cov_ctx))){
      # approach1 -------------------------------------------------------------------
      ## create empty lists for outputs 
      smooth_plot_sex_combined <- smooth_plot_sex_stratified <- vector(mode="list",length=length(roi_ctx_measurements))
      names(smooth_plot_sex_combined) <- names(smooth_plot_sex_stratified) <- roi_ctx_measurements
      assoc_res_adj = c()
      
      ## run fitting base models
      glob_ylab <- paste('adjusted brain phenotype')
      glob_xlab <- paste("adjusted BMI")
      for(roii in roi_ctx_measurements){
        analdati = subset(analdat, select=c("IID","E4_status","BMI.adj","sex",paste(roii,"adj",sep=".")))
        analdati[['y']] <- analdati[[paste(roii,"adj",sep=".")]]
        analdati = analdati %>% filter(!is.na(E4_status))
        if(do.smooth.plot){
          # APOE4-specific BMI  effect
          smooth_plot_sex_combined[[roii]] <- 
            ggplot(analdati,aes(x=BMI.adj,y=y,color=E4_status,fill=E4_status)) + 
            geom_smooth(alpha=0.2,method = 'gam', formula = y~s(x,bs='cs',k=4)) + 
            # geom_point(alpha=0.1) + 
            ggtitle(unlist(str_split(roii,"_"))[2]) +
            + xlab("")
          
          smooth_plot_sex_stratified[[roii]] <- analdati %>% 
            ggplot(aes(x=BMI.adj,y=y,color=E4_status,fill=E4_status,linetype=sex)) + 
            geom_smooth(alpha=0.2,method = 'gam', formula = y~s(x,bs='cs',k=4)) + 
            # geom_point(alpha=0.1) + 
            ggtitle(unlist(str_split(roii,"_"))[2]) +
            scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) +
            scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) + xlab("")
          
          if(str_detect(roii,"area")){
            smooth_plot_sex_stratified[[roii]] = smooth_plot_sex_stratified[[roii]] + xlab(glob_xlab) 
          }
        }
        
        if(length(unique(analdati$sex))==2){
          mod.E4.status = "y ~ BMI.adj+sex"
        }else{
          mod.E4.status = "y ~ BMI.adj"
        }
        
        ## E4-carrier
        fit.assoc.E4.carrier = lm(as.formula(mod.E4.status),data=analdati,subset=E4_status=="E4-carrier", 
                                  na.action = na.exclude)
        excl_col = which(colnames(summary(fit.assoc.E4.carrier)$coef)=="t value")
        assoc_res_adj_roii.E4.carrier = summary(fit.assoc.E4.carrier)$coef["BMI.adj",-excl_col,drop=F]
        terms = rownames(assoc_res_adj_roii.E4.carrier)
        rownames(assoc_res_adj_roii.E4.carrier)=NULL
        assoc_res_adj_roii.E4.carrier = data.frame(term=terms,assoc_res_adj_roii.E4.carrier)
        assoc_res_adj_roii.E4.carrier$N = nobs(fit.assoc.E4.carrier)
        assoc_res_adj_roii.E4.carrier$E4_status ="E4-carrier"
        
        ## E4-noncarrier
        fit.assoc.E4.noncarrier = lm(as.formula(mod.E4.status),data=analdati,subset=E4_status=="E4-noncarrier", 
                                     na.action = na.exclude)
        excl_col = which(colnames(summary(fit.assoc.E4.noncarrier)$coef)=="t value")
        assoc_res_adj_roii.E4.noncarrier = summary(fit.assoc.E4.noncarrier)$coef["BMI.adj",-excl_col,drop=F]
        terms = rownames(assoc_res_adj_roii.E4.noncarrier)
        rownames(assoc_res_adj_roii.E4.noncarrier)=NULL
        assoc_res_adj_roii.E4.noncarrier = data.frame(term=terms,assoc_res_adj_roii.E4.noncarrier)
        assoc_res_adj_roii.E4.noncarrier$N = nobs(fit.assoc.E4.noncarrier)
        assoc_res_adj_roii.E4.noncarrier$E4_status ="E4-noncarrier"
        
        # E4_status-by-BMI interaction effect
        if(length(unique(analdati$sex))==2){
          mod.int.E4.status = "y ~ E4_status*BMI.adj+sex"
        }else{
          mod.int.E4.status = "y ~ E4_status*BMI.adj"
        }
        
        fit.assoc = lm(as.formula(mod.int.E4.status),data=analdati,na.action = na.exclude)
        excl_col = which(colnames(summary(fit.assoc)$coef)=="t value")
        assoc_res_adj_roii = summary(fit.assoc)$coef[,-excl_col,drop=F]
        terms = rownames(assoc_res_adj_roii)
        rownames(assoc_res_adj_roii)=NULL
        assoc_res_adj_roii = data.frame(term=terms,assoc_res_adj_roii)
        assoc_res_adj_roii$N = nobs(fit.assoc)
        assoc_res_adj_roii$E4_status = 'APOE4-by-BMI interaction'
        
        # sex-combined BMI-effect
        res.roii = rbind(assoc_res_adj_roii.E4.carrier,
                         assoc_res_adj_roii.E4.noncarrier,
                         assoc_res_adj_roii)
        res.roii$roi = roii
        assoc_res_adj = rbind(assoc_res_adj,res.roii)
      }
      
      if(do.smooth.plot){
        ## plotting smooth curves roi-specific BMI-ctx_phenos
        sm_plot_dir = file.path(outdir,outdir2,paste("smooth_plots_adj",model,"_",filter_group,sep=''))
        if(file.exists(sm_plot_dir)){
          unlink(sm_plot_dir,recursive = T)
        }
        
        dir.create(sm_plot_dir)
        for (roi in rois){
          if(roi!="cuneus"){
            smooth_plot_roi = c(smooth_plot_sex_combined %>% tidyselect:::select(contains(roi)),
                                smooth_plot_sex_stratified %>% tidyselect:::select(contains(roi)))
          }else{
            smooth_plot_roi = c(smooth_plot_sex_combined %>% tidyselect:::select(starts_with(roi)),
                                smooth_plot_sex_stratified %>% tidyselect:::select(starts_with(roi)))
            
          }
          
          p_ranges_y <- c(ggplot_build(smooth_plot_roi[[1]])$layout$panel_scales_y[[1]]$range$range,
                          ggplot_build(smooth_plot_roi[[2]])$layout$panel_scales_y[[1]]$range$range,
                          ggplot_build(smooth_plot_roi[[3]])$layout$panel_scales_y[[1]]$range$range,
                          ggplot_build(smooth_plot_roi[[4]])$layout$panel_scales_y[[1]]$range$range,
                          ggplot_build(smooth_plot_roi[[5]])$layout$panel_scales_y[[1]]$range$range,
                          ggplot_build(smooth_plot_roi[[6]])$layout$panel_scales_y[[1]]$range$range)
          sm_roi = patchwork::wrap_plots(smooth_plot_roi) + 
            plot_layout(ncol=3,guides="collect") & 
            ylab(NULL) & 
            plot_annotation(title=paste(roi,"_",model,"Model_",filter_group,sep=""))
          
          sm_roi_common_ylim = sm_roi & 
            coord_cartesian(ylim=range(p_ranges_y))
          
          pfile = file.path(sm_plot_dir,paste(roi,"%03d.png",sep=""))
          png(file=pfile,width=11,height=8,units='in',res=300)
          print(sm_roi)
          print(sm_roi_common_ylim)
          dev.off()
        }#for (roi in rois)
      }
      rm(analdat)
      
      # approach2: covaraite adjustment after splitting -----------------------
      analdat = d %>% mutate(age.c = scale(age)[,1],
                             age.c2 = (scale(age)/2)^2) %>% # for stability (Gelman)
        mutate(logBMI = log(BMI), ICV.scaled = scale(ICV)[,1]/2)
      analdat = get_subset("analdat",filter_expression_i)
      
      ## removing outliers of BMI and ICV
      d_cleaned = remove.outliers_grubbs(analdat,varnames = c('logBMI',"ICV.scaled"))
      analdat = d_cleaned$cleaned_data
      print(d_cleaned$counts.NA)
      rm(d_cleaned)
      
      covariate_names = setdiff(c(cov_ctx,cov_BMI),c("age","sex"))
      assoc_res_cov = c()
      for(roii in roi_ctx_measurements){
        analdati = subset(analdat,select=c("IID",roii,'logBMI','ICV.scaled','E4_status','sex',
                                           setdiff(c(cov_ctx,cov_BMI),c("age","sex","BMI","ICV")))) 
        
        analdati[['y']] <- inormal(analdati[[paste(roii,sep=".")]])
        analdati[['BMI']] <- inormal(analdati[['logBMI']])
        analdati[['ICV']] <- inormal(analdati[['ICV.scaled']])
        covariate_names = covariate_names[which(!psych::describe(subset(analdati,select=covariate_names))$sd==0)]
        
        mod0 = paste0(c(covariate_names),collapse = " + ")
        if(length(unique(analdati$sex))==2){
          mod0 = paste("y ~ E4_status*BMI + sex*(",mod0,")", sep="")
        }else{#length(unique(analdati$sex))=1
          mod0 = paste("y ~ E4_status*BMI + ", mod0, sep="")
        }
        capture.output(mod0,file=file.path(outdir,input_specification_file), append=T)
        
        fit.assoc = lm(mod0,data=analdati,na.action = na.exclude)
        assoc_res_adj_roii = summary(fit.assoc)$coef
        terms = rownames(assoc_res_adj_roii)
        excl_col = which(colnames(assoc_res_adj_roii)=="t value")
        assoc_res_adj_roii = data.frame(term=terms,assoc_res_adj_roii[,-excl_col])
        assoc_res_adj_roii = subset(assoc_res_adj_roii, term %in% c("BMI","E4_statusE4-carrier:BMI"))
        assoc_res_adj_roii$N = nobs(fit.assoc)
        assoc_res_adj_roii$roi = roii
        assoc_res_adj_roii$E4_status = 'APOE4-by-BMI interaction'
        assoc_res_adj_roii$mod = mod0
        
        rownames(assoc_res_adj_roii)=NULL
        assoc_res_cov = rbind(assoc_res_cov,assoc_res_adj_roii)
      }
      
      # formatting results ----------------------------------------------------------#
      res1 = subset(assoc_res_adj,term %in% c("BMI.adj","E4_status:BMI.adj"))
      names(res1)[2:4] = c("Estimate","SE","P")
      res1$method = 'm1'
      res1$term = str_remove(res1$term,"[.]adj")
      
      res2 = assoc_res_cov
      names(res2)[2:4] = c("Estimate","SE","P")
      res2$method = 'm2'
      
      # creating outputs ------------------------------------------------------------#
      # table containing both sets of the results
      roi_assoc_res_two_methods = merge(subset(res1,E4_status=="APOE4-by-BMI interaction"),
                                        res2, by=c('roi',"term"))
      names(roi_assoc_res_two_methods) = str_replace(names(roi_assoc_res_two_methods),"[.]x","_reg_out_covs")
      names(roi_assoc_res_two_methods) = str_replace(names(roi_assoc_res_two_methods),"[.]y","_adj_covs_in_mod")
      
      # cleaning-up -----------------------------------------------------------------#
      save(
        res1,res2,roi_assoc_res_two_methods,model,group_name,filter_expression_i,filter_group,
        file=file.path(outdir,outdir2,paste("roi_assoc_res_two_methods_",model,"_",filter_group,".Rdata",sep=''))
      )
      rm(res1,res2,res_two_methods,adj_ctxTH,adj_BMI,analdat,analdati,group_name,filter_expression_i)
    }#if(min(table(analdat$E4_status)) >(30+length(cov_ctx)))
    
    #ending message ---------------------------------------------------------------#
    end.msg = paste("\nFinished obtaining association estimates for BMI vs. roi-ctx phenotypes under ",model," model (",filter_group,").\n",sep='')
    cat(end.msg)
    rm(filter_group)
    sink()
    closeAllConnections()
    print(readLines(logfile))
  } #for(filter_i)
  rm(model)
} #for(model in models)
options(op)

## ----forest_plots_groups----------------------------------------------------------------------------
bmi.strata.levels=
  # all
  c('BMI <20kg/m2',
    "20<= BMI <25kg/m2",
    "25<= BMI <30kg/m2", 
    "BMI >=30kg/m2",
    # youger
    "age <65 years & BMI <20kg/m2", 
    "age <65 years & 20<= BMI <25kg/m2", 
    "age <65 years & 25<= BMI <30kg/m2", 
    "age <65 years & BMI >=30kg/m2", 
    # older
    "age >=65 years & BMI <20kg/m2", 
    "age >=65 years & 20<= BMI <25kg/m2", 
    "age >=65 years & 25<= BMI <30kg/m2", 
    "age >=65 years & BMI >=30kg/m2")


## ----forest_plot_all--------------------------------------------------------------------------------
# in all participants
sel_bmi_levels=c(2:4); age_group="all"
filter_is=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$index 

# sample size table
a=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$group_name
b=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$N_E4_noncarrier
c=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$N_E4_carrier
sample.sizes = data.table(rbind(b,c))
names(sample.sizes) = a
sample.sizes

# format association results to be plotted (e.g., add 95CI)
res1_all=c()
ranges_Estimate=c()
ranges_95CI=c()
for (model_i in c('base')){
  for(filter_i in filter_is){
    sampleN_E4noncarrier=df_filter$N_E4_noncarrier[filter_i]
    sampleN_E4carrier=df_filter$N_E4_carrier[filter_i]
    group_name=df_filter$group_name[filter_i]
    filter_group=df_filter$group_label[filter_i]
    filter_expression_i=df_filter$expression[filter_i]
    if(is.null(outdir2)){
      load(file.path(outdir,paste("roi_assoc_res_two_methods_",model_i,"_",filter_group,".Rdata",sep='')))
    }else{
      load(file.path(outdir,outdir2,paste("roi_assoc_res_two_methods_",model_i,"_",filter_group,".Rdata",sep='')))
    }
    res1_all=rbind(res1_all,res1 %>% mutate(model=model_i, group=group_name, 
                                            Nmax_E4noncarrier=sampleN_E4noncarrier,
                                            Nmax_E4carrier=sampleN_E4carrier,
                                            U95CI=Estimate + 1.96*SE,L95CI=Estimate - 1.96*SE))
    
    res1_all = res1_all %>% filter(str_detect(roi,"_thickness"))
    
    ranges_Estimate=rbind(ranges_Estimate, range(res1_all$Estimate))
    ranges_95CI=rbind(ranges_95CI, range(res1_all$L95CI,res1_all$U95CI))
  }
}

res1_all=res1_all %>% 
  mutate(group=factor(group, 
                      levels=bmi.strata.levels))
data.frame(table(res1_all$group))

plotd=res1_all %>% 
  # mutate(U95CI=Estimate + 1.96*SE,L95CI=Estimate - 1.96*SE) %>%
  mutate(pheno=str_split(res1_all$roi,"_",simplify=T)[,2]) %>%
  mutate(roi=str_split(res1_all$roi,"_",simplify=T)[,1]) %>%
  filter(E4_status != 'APOE4-by-BMI interaction') %>% 
  mutate(E4_status=factor(E4_status,levels=c("E4-noncarrier","E4-carrier"))) %>%
  mutate(group=factor(group, levels=(bmi.strata.levels)))
yrange_forest=0.85*c(-1,1)*max(abs(ranges_95CI))

cleanup=theme(panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              panel.grid.major.y=element_line(colour="grey90",linetype=3),
              panel.background=element_blank(),
              axis.line=element_line(colour="black"),
              panel.border=element_rect(colour="grey", fill=NA, linewidth=0.5))

## sample sizes ---------------------------------------------------------------#
plotd %>% ggplot(aes(x=group, y=N, color=E4_status, shape=model)) + geom_point() + 
  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))

## plotting -------------------------------------------------------------------#
plot.alpha=0.7
pos <- position_dodge(width=0)#

# for (model_plot in c("base","full")){
for (model_plot in c("base")){
  for (pheno_plot in c("thickness") ){
    plotd_i=plotd %>% filter(model==model_plot, pheno==pheno_plot) %>% arrange(roi,group)
    roi.levels=plotd_i %>% filter(E4_status=="E4-carrier",group==sort(unique(plotd_i$group))[length(sort(unique(plotd_i$group)))]) %>% arrange(Estimate)
    roi.levels=roi.levels$roi
    #################################################
    roi.levels = fread('roi_levels.txt',header=F)$V1
    roi.levels = rev(roi.levels)
    #################################################
    plotd_i=plotd_i %>% mutate(roi=factor(roi,levels=rev(roi.levels))) %>% arrange(roi)
    
    ## 1
    p <- plotd_i %>%
      ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,
                 group=group,color=group,fill=group,
                 shape=group)) + 
      geom_point(size=1.15, stroke=1, position=pos,alpha=plot.alpha) +
      scale_shape_manual(values=c(1,2,6,12)) + 
      geom_hline(yintercept=0, size=0.4) +
      geom_path(position=pos,alpha=0.5) +
      geom_hline(yintercept=0, linewidth=0.5)+
      xlab(NULL) + 
      coord_flip(ylim=yrange_forest) +
      scale_y_continuous(labels=scales::label_number(accuracy=0.01)) + 
      ggtitle(paste("BMI vs. cortical",pheno_plot,"associations")) + 
      # shape= guide_legend(reverse=F,title=element_blank()),
      # fill=guide_legend(show=FALSE) ) +
      geom_ribbon(aes(ymin=L95CI,ymax=U95CI), alpha=0.15, position=pos, color=NA) #+       guides(col=guide_legend(reverse=F,title=element_blank()), 
    
    p = p+
      facet_grid(.~E4_status)+
      scale_fill_manual(values = rev(c("#e7298a","#7570b3","#d95f02")))+
      scale_color_manual(values = rev(c("#e7298a","#7570b3","#d95f02")))
    # scale_fill_brewer(palette="Dark2")+
    # scale_color_brewer(palette="Dark2") 
    p = p+
      cleanup+ 
      labs(color=NULL,fill=NULL)+
      labs(shape=NULL)+
      theme(legend.position='top',legend.box="vertical", legend.margin=margin())+ 
      # theme(plot.margin=margin(t=0.15, r=1.5, b=0.15, l=0.15, "cm"),
      #       axis.text.x=element_text(angle=-30, hjust=0, vjust=0.5)) + 
      ylab("Estimate in SD units")
    print(p)
    
    pfile=paste("forest_plot_by_E4_status_",
                age_group,"_",
                model_plot,"_",
                pheno_plot,"_",Sys.Date(),".png",sep='')
    if(is.null(outdir2)){
      pfile=file.path(outdir,pfile)
    }else{
      pfile=file.path(outdir,outdir2,pfile)
    }
    png(pfile,width=12*0.8,height=8.5*0.8,units='in',res=300)
    print(p)
    dev.off()
    
    ## by bmi classes
    p_BMI <- plotd_i %>%
      ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,
                 group=E4_status,color=E4_status,fill=E4_status,
                 shape=group)) + 
      geom_point(size=1.15, stroke=1, position=pos,alpha=plot.alpha) +
      scale_shape_manual(values=c(1,2,6,12)) + 
      geom_hline(yintercept=0, size=0.4) +
      geom_path(position=pos,alpha=0.5) +
      geom_hline(yintercept=0, linewidth=0.5)+
      xlab(NULL) + 
      coord_flip(ylim=yrange_forest) +
      scale_y_continuous(labels=scales::label_number(accuracy=0.01)) + 
      ggtitle(paste("BMI vs. cortical",pheno_plot,"associations")) + 
      geom_ribbon(aes(ymin=L95CI,ymax=U95CI), alpha=0.15, position=pos, color=NA) #+       guides(col=guide_legend(reverse=F,title=element_blank()), 
    # shape= guide_legend(reverse=F,title=element_blank()),
    # fill=guide_legend(show=FALSE) ) +
    
    # by E4-category
    p_BMI = p_BMI+ 
      facet_grid(.~group)+scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
      scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))
    p_BMI+
      cleanup+ 
      labs(color=NULL,fill=NULL)+
      labs(shape=NULL)+
      theme(legend.position='top',legend.box="vertical", legend.margin=margin())+ 
      # theme(plot.margin=margin(t=0.15, r=1.5, b=0.15, l=0.15, "cm"),
      #       axis.text.x=element_text(angle=-30, hjust=0, vjust=0.5)) + 
      ylab("Estimate in SD units")
    
    pfile_BMI=paste("forest_plot_by_BMI_classes_",
                    age_group,"_",
                    model_plot,"_",
                    pheno_plot,"_",Sys.Date(),".png",sep='')
    if(is.null(outdir2)){
      pfile_BMI=file.path(outdir,pfile_BMI)
    }else{
      pfile_BMI=file.path(outdir,outdir2,pfile_BMI)
    }
    ggsave(pfile_BMI,width=12*0.8,height=8.5*0.8,units='in',dpi=300)
    
  }
}


## ----forest_plot_younger----------------------------------------------------------------------------
# in younger participants
sel_bmi_levels=c(6:8); age_group="age_lt_65"

filter_is=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$index 
# sample size table
a=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$group_name
b=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$N_E4_noncarrier
c=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$N_E4_carrier
sample.sizes = data.table(rbind(b,c))
names(sample.sizes) = a
print(sample.sizes)

# format association results to be plotted (e.g., add 95CI)
res1_all=c()
ranges_Estimate=c()
ranges_95CI=c()
for (model_i in c('base')){
  for(filter_i in filter_is){
    sampleN_E4noncarrier=df_filter$N_E4_noncarrier[filter_i]
    sampleN_E4carrier=df_filter$N_E4_carrier[filter_i]
    group_name=df_filter$group_name[filter_i]
    filter_group=df_filter$group_label[filter_i]
    filter_expression_i=df_filter$expression[filter_i]
    if(is.null(outdir2)){
      load(file.path(outdir,paste("roi_assoc_res_two_methods_",model_i,"_",filter_group,".Rdata",sep='')))
    }else{
      load(file.path(outdir,outdir2,paste("roi_assoc_res_two_methods_",model_i,"_",filter_group,".Rdata",sep='')))
    }
    res1_all=rbind(res1_all,res1 %>% mutate(model=model_i, group=group_name, 
                                            Nmax_E4noncarrier=sampleN_E4noncarrier,
                                            Nmax_E4carrier=sampleN_E4carrier,
                                            U95CI=Estimate + 1.96*SE,L95CI=Estimate - 1.96*SE))
    
    res1_all = res1_all %>% filter(str_detect(roi,"_thickness"))
    
    ranges_Estimate=rbind(ranges_Estimate, range(res1_all$Estimate))
    ranges_95CI=rbind(ranges_95CI, range(res1_all$L95CI,res1_all$U95CI))
  }
}

res1_all=res1_all %>% 
  mutate(group=factor(group, 
                      levels=bmi.strata.levels))
data.frame(table(res1_all$group))

plotd=res1_all %>% 
  # mutate(U95CI=Estimate + 1.96*SE,L95CI=Estimate - 1.96*SE) %>%
  mutate(pheno=str_split(res1_all$roi,"_",simplify=T)[,2]) %>%
  mutate(roi=str_split(res1_all$roi,"_",simplify=T)[,1]) %>%
  filter(E4_status != 'APOE4-by-BMI interaction') %>% 
  mutate(E4_status=factor(E4_status,levels=c("E4-noncarrier","E4-carrier"))) %>%
  mutate(group=factor(group, levels=(bmi.strata.levels)))
yrange_forest=0.85*c(-1,1)*max(abs(ranges_95CI))

cleanup=theme(panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              panel.grid.major.y=element_line(colour="grey90",linetype=3),
              panel.background=element_blank(),
              axis.line=element_line(colour="black"),
              panel.border=element_rect(colour="grey", fill=NA, linewidth=0.5))

## sample sizes ---------------------------------------------------------------#
plotd %>% ggplot(aes(x=group, y=N, color=E4_status, shape=model)) + geom_point() + 
  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))

## plotting -------------------------------------------------------------------#
plot.alpha=0.7
pos <- position_dodge(width=0)#

# for (model_plot in c("base","full")){
for (model_plot in c("base")){
  for (pheno_plot in c("thickness") ){
    plotd_i=plotd %>% filter(model==model_plot, pheno==pheno_plot) %>% arrange(roi,group)
    roi.levels=plotd_i %>% filter(E4_status=="E4-carrier",group==sort(unique(plotd_i$group))[length(sort(unique(plotd_i$group)))]) %>% arrange(Estimate)
    roi.levels=roi.levels$roi
    plotd_i=plotd_i %>% mutate(roi=factor(roi,levels=rev(roi.levels))) %>% arrange(roi)
    
    ## 1
    p <- plotd_i %>%
      ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,
                 group=group,color=group,fill=group,
                 shape=group)) + 
      geom_point(size=1.15, stroke=1, position=pos,alpha=plot.alpha) +
      scale_shape_manual(values=c(1,2,6,12)) + 
      geom_hline(yintercept=0, size=0.4) +
      geom_path(position=pos,alpha=0.5) +
      geom_hline(yintercept=0, linewidth=0.5)+
      xlab(NULL) + 
      coord_flip(ylim=yrange_forest) +
      scale_y_continuous(labels=scales::label_number(accuracy=0.01)) + 
      ggtitle(paste("BMI vs. cortical",pheno_plot,"associations")) + 
      # shape= guide_legend(reverse=F,title=element_blank()),
      # fill=guide_legend(show=FALSE) ) +
      geom_ribbon(aes(ymin=L95CI,ymax=U95CI), alpha=0.15, position=pos, color=NA) #+       guides(col=guide_legend(reverse=F,title=element_blank()), 
    
    p = p+
      facet_grid(.~E4_status)+
      scale_fill_manual(values = rev(c("#e7298a","#7570b3","#d95f02")))+
      scale_color_manual(values = rev(c("#e7298a","#7570b3","#d95f02")))
    # scale_fill_brewer(palette="Dark2")+
    # scale_color_brewer(palette="Dark2") 
    p+
      cleanup+ 
      labs(color=NULL,fill=NULL)+
      labs(shape=NULL)+
      theme(legend.position='top',legend.box="vertical", legend.margin=margin())+ 
      # theme(plot.margin=margin(t=0.15, r=1.5, b=0.15, l=0.15, "cm"),
      #       axis.text.x=element_text(angle=-30, hjust=0, vjust=0.5)) + 
      ylab("Estimate in SD units")
    
    pfile=paste("forest_plot_by_E4_status_",
                age_group,"_",
                model_plot,"_",
                pheno_plot,"_",Sys.Date(),".png",sep='')
    if(is.null(outdir2)){
      pfile=file.path(outdir,pfile)
    }else{
      pfile=file.path(outdir,outdir2,pfile)
    }
    ggsave(pfile,width=12*0.8,height=8.5*0.8,units='in',dpi=300)
    
    ## by bmi classes
    p_BMI <- plotd_i %>%
      ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,
                 group=E4_status,color=E4_status,fill=E4_status,
                 shape=group)) + 
      geom_point(size=1.15, stroke=1, position=pos,alpha=plot.alpha) +
      scale_shape_manual(values=c(1,2,6,12)) + 
      geom_hline(yintercept=0, size=0.4) +
      geom_path(position=pos,alpha=0.5) +
      geom_hline(yintercept=0, linewidth=0.5)+
      xlab(NULL) + 
      coord_flip(ylim=yrange_forest) +
      scale_y_continuous(labels=scales::label_number(accuracy=0.01)) + 
      ggtitle(paste("BMI vs. cortical",pheno_plot,"associations")) + 
      geom_ribbon(aes(ymin=L95CI,ymax=U95CI), alpha=0.15, position=pos, color=NA) #+       guides(col=guide_legend(reverse=F,title=element_blank()), 
    # shape= guide_legend(reverse=F,title=element_blank()),
    # fill=guide_legend(show=FALSE) ) +
    
    # by E4-category
    p_BMI = p_BMI+ 
      facet_grid(.~group)+scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
      scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))
    p_BMI = p_BMI+
      cleanup+ 
      labs(color=NULL,fill=NULL)+
      labs(shape=NULL)+
      theme(legend.position='top',legend.box="vertical", legend.margin=margin())+ 
      # theme(plot.margin=margin(t=0.15, r=1.5, b=0.15, l=0.15, "cm"),
      #       axis.text.x=element_text(angle=-30, hjust=0, vjust=0.5)) + 
      ylab("Estimate in SD units")
    
    pfile_BMI=paste("forest_plot_by_BMI_classes_",
                    age_group,"_",
                    model_plot,"_",
                    pheno_plot,"_",Sys.Date(),".png",sep='')
    if(is.null(outdir2)){
      pfile_BMI=file.path(outdir,pfile_BMI)
    }else{
      pfile_BMI=file.path(outdir,outdir2,pfile_BMI)
    }
    print(p_BMI)
    
    png(pfile_BMI,width=12*0.8,height=8.5*0.8,units='in',res=300)
    print(p_BMI)
    dev.off()
    
  }
}


## ----forest_plot_older-----------------------------------------------------------------------------
# in older participants
sel_bmi_levels=c(10:12); age_group="age_gt_65"

filter_is=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$index 
# sample size table
a=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$group_name
b=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$N_E4_noncarrier
c=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$N_E4_carrier
sample.sizes = data.table(rbind(b,c))
names(sample.sizes) = a
sample.sizes

# format association results to be plotted (e.g., add 95CI)
res1_all=c()
ranges_Estimate=c()
ranges_95CI=c()
for (model_i in c('base')){
  for(filter_i in filter_is){
    sampleN_E4noncarrier=df_filter$N_E4_noncarrier[filter_i]
    sampleN_E4carrier=df_filter$N_E4_carrier[filter_i]
    group_name=df_filter$group_name[filter_i]
    filter_group=df_filter$group_label[filter_i]
    filter_expression_i=df_filter$expression[filter_i]
    if(is.null(outdir2)){
      load(file.path(outdir,paste("roi_assoc_res_two_methods_",model_i,"_",filter_group,".Rdata",sep='')))
    }else{
      load(file.path(outdir,outdir2,paste("roi_assoc_res_two_methods_",model_i,"_",filter_group,".Rdata",sep='')))
    }
    res1_all=rbind(res1_all,res1 %>% mutate(model=model_i, group=group_name, 
                                            Nmax_E4noncarrier=sampleN_E4noncarrier,
                                            Nmax_E4carrier=sampleN_E4carrier,
                                            U95CI=Estimate + 1.96*SE,L95CI=Estimate - 1.96*SE))
    ranges_Estimate=rbind(ranges_Estimate, range(res1_all$Estimate))
    
    res1_all = res1_all %>% filter(str_detect(roi,"_thickness"))
    
    ranges_95CI=rbind(ranges_95CI, range(res1_all$L95CI,res1_all$U95CI))
  }
}

res1_all=res1_all %>% 
  mutate(group=factor(group, 
                      levels=bmi.strata.levels))
data.frame(table(res1_all$group))

plotd=res1_all %>% 
  # mutate(U95CI=Estimate + 1.96*SE,L95CI=Estimate - 1.96*SE) %>%
  mutate(pheno=str_split(res1_all$roi,"_",simplify=T)[,2]) %>%
  mutate(roi=str_split(res1_all$roi,"_",simplify=T)[,1]) %>%
  filter(E4_status != 'APOE4-by-BMI interaction') %>% 
  mutate(E4_status=factor(E4_status,levels=c("E4-noncarrier","E4-carrier"))) %>%
  mutate(group=factor(group, levels=(bmi.strata.levels)))
yrange_forest=0.75*c(-1,1)*max(abs(ranges_95CI))

cleanup=theme(panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              panel.grid.major.y=element_line(colour="grey90",linetype=3),
              panel.background=element_blank(),
              axis.line=element_line(colour="black"),
              panel.border=element_rect(colour="grey", fill=NA, linewidth=0.5))

## sample sizes ---------------------------------------------------------------#
plotd %>% ggplot(aes(x=group, y=N, color=E4_status, shape=model)) + geom_point() + 
  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))

## plotting -------------------------------------------------------------------#
plot.alpha=0.7
pos <- position_dodge(width=0)#

# for (model_plot in c("base","full")){
for (model_plot in c("base")){
  for (pheno_plot in c("thickness") ){
    plotd_i=plotd %>% filter(model==model_plot, pheno==pheno_plot) %>% arrange(roi,group)
    roi.levels=plotd_i %>% filter(E4_status=="E4-carrier",group==sort(unique(plotd_i$group))[length(sort(unique(plotd_i$group)))]) %>% arrange(Estimate)
    roi.levels=roi.levels$roi
    plotd_i=plotd_i %>% mutate(roi=factor(roi,levels=rev(roi.levels))) %>% arrange(roi)
    
    ## 1
    p <- plotd_i %>%
      ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,
                 group=group,color=group,fill=group,
                 shape=group)) + 
      geom_point(size=1.15, stroke=1, position=pos,alpha=plot.alpha) +
      scale_shape_manual(values=c(1,2,6,12)) + 
      geom_hline(yintercept=0, size=0.4) +
      geom_path(position=pos,alpha=0.5) +
      geom_hline(yintercept=0, linewidth=0.5)+
      xlab(NULL) + 
      coord_flip(ylim=yrange_forest) +
      scale_y_continuous(labels=scales::label_number(accuracy=0.01)) + 
      ggtitle(paste("BMI vs. cortical",pheno_plot,"associations")) + 
      # shape= guide_legend(reverse=F,title=element_blank()),
      # fill=guide_legend(show=FALSE) ) +
      geom_ribbon(aes(ymin=L95CI,ymax=U95CI), alpha=0.15, position=pos, color=NA) #+       guides(col=guide_legend(reverse=F,title=element_blank()), 
    
    p = p+
      facet_grid(.~E4_status)+
      scale_fill_manual(values = rev(c("#e7298a","#7570b3","#d95f02")))+
      scale_color_manual(values = rev(c("#e7298a","#7570b3","#d95f02")))
    # scale_fill_brewer(palette="Dark2")+
    # scale_color_brewer(palette="Dark2") 
    p+
      cleanup+ 
      labs(color=NULL,fill=NULL)+
      labs(shape=NULL)+
      theme(legend.position='top',legend.box="vertical", legend.margin=margin())+ 
      # theme(plot.margin=margin(t=0.15, r=1.5, b=0.15, l=0.15, "cm"),
      #       axis.text.x=element_text(angle=-30, hjust=0, vjust=0.5)) + 
      ylab("Estimate in SD units")
    
    pfile=paste("forest_plot_by_E4_status_",
                age_group,"_",
                model_plot,"_",
                pheno_plot,"_",Sys.Date(),".png",sep='')
    if(is.null(outdir2)){
      pfile=file.path(outdir,pfile)
    }else{
      pfile=file.path(outdir,outdir2,pfile)
    }
    ggsave(pfile,width=12*0.8,height=8.5*0.8,units='in',dpi=300)
    
    ## by bmi classes
    p_BMI <- plotd_i %>%
      ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,
                 group=E4_status,color=E4_status,fill=E4_status,
                 shape=group)) + 
      geom_point(size=1.15, stroke=1, position=pos,alpha=plot.alpha) +
      scale_shape_manual(values=c(1,2,6,12)) + 
      geom_hline(yintercept=0, size=0.4) +
      geom_path(position=pos,alpha=0.5) +
      geom_hline(yintercept=0, linewidth=0.5)+
      xlab(NULL) + 
      coord_flip(ylim=yrange_forest) +
      scale_y_continuous(labels=scales::label_number(accuracy=0.01)) + 
      ggtitle(paste("BMI vs. cortical",pheno_plot,"associations")) + 
      geom_ribbon(aes(ymin=L95CI,ymax=U95CI), alpha=0.15, position=pos, color=NA) #+       guides(col=guide_legend(reverse=F,title=element_blank()), 
    # shape= guide_legend(reverse=F,title=element_blank()),
    # fill=guide_legend(show=FALSE) ) +
    
    # by E4-category
    p_BMI = p_BMI+ 
      facet_grid(.~group)+scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
      scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))
    p_BMI = p_BMI+
      cleanup+ 
      labs(color=NULL,fill=NULL)+
      labs(shape=NULL)+
      theme(legend.position='top',legend.box="vertical", legend.margin=margin())+ 
      # theme(plot.margin=margin(t=0.15, r=1.5, b=0.15, l=0.15, "cm"),
      #       axis.text.x=element_text(angle=-30, hjust=0, vjust=0.5)) + 
      ylab("Estimate in SD units")
    
    pfile_BMI=paste("forest_plot_by_BMI_classes_",
                    age_group,"_",
                    model_plot,"_",
                    pheno_plot,"_",Sys.Date(),".png",sep='')
    if(is.null(outdir2)){
      pfile_BMI=file.path(outdir,pfile_BMI)
    }else{
      pfile_BMI=file.path(outdir,outdir2,pfile_BMI)
    }
    print(p_BMI)
    png(pfile_BMI,width=12*0.8,height=8.5*0.8,units='in',res=300)
    print(p_BMI)
    dev.off()
  }
}


## ----forest_plot_bmi-pooled-------------------------------------------------------------------------
# in bmi-pooled participants
bmi.strata.levels=
  # bmi-pooled
  c('all',
    "age <65 years",
    "age >=65 years")

sel_bmi_levels=c(1,2,3); age_group="bmi-pooled"
filter_is=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$index 

# sample size table
a=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$group_name
b=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$N_E4_noncarrier
c=(df_filter %>% filter(group_name %in% bmi.strata.levels[sel_bmi_levels]))$N_E4_carrier
sample.sizes = data.table(rbind(b,c))
names(sample.sizes) = a
sample.sizes

# format association results to be plotted (e.g., add 95CI)
res1_all=c()
ranges_Estimate=c()
ranges_95CI=c()
for (model_i in c('base')){
  for(filter_i in filter_is){
    sampleN_E4noncarrier=df_filter$N_E4_noncarrier[filter_i]
    sampleN_E4carrier=df_filter$N_E4_carrier[filter_i]
    group_name=df_filter$group_name[filter_i]
    filter_group=df_filter$group_label[filter_i]
    filter_expression_i=df_filter$expression[filter_i]
    if(is.null(outdir2)){
      load(file.path(outdir,paste("roi_assoc_res_two_methods_",model_i,"_",filter_group,".Rdata",sep='')))
    }else{
      load(file.path(outdir,outdir2,paste("roi_assoc_res_two_methods_",model_i,"_",filter_group,".Rdata",sep='')))
    }
    res1_all=rbind(res1_all,res1 %>% mutate(model=model_i, group=group_name, 
                                            Nmax_E4noncarrier=sampleN_E4noncarrier,
                                            Nmax_E4carrier=sampleN_E4carrier,
                                            U95CI=Estimate + 1.96*SE,L95CI=Estimate - 1.96*SE))
    
    res1_all = res1_all %>% filter(str_detect(roi,"_thickness"))
    
    ranges_Estimate=rbind(ranges_Estimate, range(res1_all$Estimate))
    ranges_95CI=rbind(ranges_95CI, range(res1_all$L95CI,res1_all$U95CI))
  }
}

res1_all=res1_all %>% 
  mutate(group=factor(group, 
                      levels=bmi.strata.levels))
data.frame(table(res1_all$group))

plotd=res1_all %>% 
  # mutate(U95CI=Estimate + 1.96*SE,L95CI=Estimate - 1.96*SE) %>%
  mutate(pheno=str_split(res1_all$roi,"_",simplify=T)[,2]) %>%
  mutate(roi=str_split(res1_all$roi,"_",simplify=T)[,1]) %>%
  filter(E4_status != 'APOE4-by-BMI interaction') %>% 
  mutate(E4_status=factor(E4_status,levels=c("E4-noncarrier","E4-carrier"))) %>%
  mutate(group=factor(group, levels=(bmi.strata.levels)))
yrange_forest=0.85*c(-1,1)*max(abs(ranges_95CI))

cleanup=theme(panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              panel.grid.major.y=element_line(colour="grey90",linetype=3),
              panel.background=element_blank(),
              axis.line=element_line(colour="black"),
              panel.border=element_rect(colour="grey", fill=NA, linewidth=0.5))

## sample sizes ---------------------------------------------------------------#
plotd %>% ggplot(aes(x=group, y=N, color=E4_status, shape=model)) + geom_point() + 
  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))

## plotting -------------------------------------------------------------------#
plot.alpha=0.7
pos <- position_dodge(width=0)#

# for (model_plot in c("base","full")){
for (model_plot in c("base")){
  for (pheno_plot in c("thickness") ){
    plotd_i=plotd %>% filter(model==model_plot, pheno==pheno_plot) %>% arrange(roi,group)
    roi.levels=plotd_i %>% filter(E4_status=="E4-carrier",group==sort(unique(plotd_i$group))[length(sort(unique(plotd_i$group)))]) %>% arrange(Estimate)
    roi.levels=roi.levels$roi
    plotd_i=plotd_i %>% mutate(roi=factor(roi,levels=rev(roi.levels))) %>% arrange(roi)
    
    ## 1
    p <- plotd_i %>%
      ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,
                 group=group,color=group,fill=group,
                 shape=group)) + 
      geom_point(size=1.15, stroke=1, position=pos,alpha=plot.alpha) +
      scale_shape_manual(values=c(1,2,6,12)) + 
      geom_hline(yintercept=0, size=0.4) +
      geom_path(position=pos,alpha=0.5) +
      geom_hline(yintercept=0, linewidth=0.5)+
      xlab(NULL) + 
      coord_flip(ylim=yrange_forest) +
      scale_y_continuous(labels=scales::label_number(accuracy=0.01)) + 
      ggtitle(paste("BMI vs. cortical",pheno_plot,"associations")) + 
      # shape= guide_legend(reverse=F,title=element_blank()),
      # fill=guide_legend(show=FALSE) ) +
      geom_ribbon(aes(ymin=L95CI,ymax=U95CI), alpha=0.15, position=pos, color=NA) #+       guides(col=guide_legend(reverse=F,title=element_blank()), 
    
    p = p+
      facet_grid(.~E4_status)+
      scale_fill_manual(values = rev(c("#e7298a","#7570b3","#d95f02")))+
      scale_color_manual(values = rev(c("#e7298a","#7570b3","#d95f02")))
    # scale_fill_brewer(palette="Dark2")+
    # scale_color_brewer(palette="Dark2") 
    p = p+
      cleanup+ 
      labs(color=NULL,fill=NULL)+
      labs(shape=NULL)+
      theme(legend.position='top',legend.box="vertical", legend.margin=margin())+ 
      # theme(plot.margin=margin(t=0.15, r=1.5, b=0.15, l=0.15, "cm"),
      #       axis.text.x=element_text(angle=-30, hjust=0, vjust=0.5)) + 
      ylab("Estimate in SD units")
    print(p)
    
    pfile=paste("forest_plot_by_E4_status_",
                age_group,"_",
                model_plot,"_",
                pheno_plot,"_",Sys.Date(),".png",sep='')
    if(is.null(outdir2)){
      pfile=file.path(outdir,pfile)
    }else{
      pfile=file.path(outdir,outdir2,pfile)
    }
    png(pfile,width=12*0.8,height=8.5*0.8,units='in',res=300)
    print(p)
    dev.off()
    
    ## by bmi classes
    p_BMI <- plotd_i %>%
      ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,
                 group=E4_status,color=E4_status,fill=E4_status,
                 shape=group)) + 
      geom_point(size=1.15, stroke=1, position=pos,alpha=plot.alpha) +
      scale_shape_manual(values=c(1,2,6,12)) + 
      geom_hline(yintercept=0, size=0.4) +
      geom_path(position=pos,alpha=0.5) +
      geom_hline(yintercept=0, linewidth=0.5)+
      xlab(NULL) + 
      coord_flip(ylim=yrange_forest) +
      scale_y_continuous(labels=scales::label_number(accuracy=0.01)) + 
      ggtitle(paste("BMI vs. cortical",pheno_plot,"associations")) + 
      geom_ribbon(aes(ymin=L95CI,ymax=U95CI), alpha=0.15, position=pos, color=NA) #+       guides(col=guide_legend(reverse=F,title=element_blank()), 
    # shape= guide_legend(reverse=F,title=element_blank()),
    # fill=guide_legend(show=FALSE) ) +
    
    # by E4-category
    p_BMI = p_BMI+ 
      facet_grid(.~group)+scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
      scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))
    p_BMI+
      cleanup+ 
      labs(color=NULL,fill=NULL)+
      labs(shape=NULL)+
      theme(legend.position='top',legend.box="vertical", legend.margin=margin())+ 
      # theme(plot.margin=margin(t=0.15, r=1.5, b=0.15, l=0.15, "cm"),
      #       axis.text.x=element_text(angle=-30, hjust=0, vjust=0.5)) + 
      ylab("Estimate in SD units")
    
    pfile_BMI=paste("forest_plot_by_BMI_classes_",
                    age_group,"_",
                    model_plot,"_",
                    pheno_plot,"_",Sys.Date(),".png",sep='')
    if(is.null(outdir2)){
      pfile_BMI=file.path(outdir,pfile_BMI)
    }else{
      pfile_BMI=file.path(outdir,outdir2,pfile_BMI)
    }
    ggsave(pfile_BMI,width=12*0.8,height=8.5*0.8,units='in',dpi=300)
    
  }
}
