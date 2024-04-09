#*****************************************************************************#
#
# 3:
# Examine associations between BMI vs. ctx-phenotype across the 34 regions
# adjusting for age, sex, ICV (for brain phenotype) and any 
# cohort-specific covariates in the participants according to the filter grouping.
#
#*****************************************************************************#
get_adj_ctx = function(roii_pheno,d,covariate_names){
  d = d %>% mutate(age.c = scale(age)[,1],
                   age.c2 = (scale(age)/2)^2)# for stability (Gelman)
  analdati = subset(d,select=c(covariate_names,'age.c', 'age.c2',roii_pheno)) 
  
  ## remove outliers in the original variable
  y_outliers_removed0 = remove.outliers_grubbs(analdati,varnames=roii_pheno)
  print(y_outliers_removed0$counts.NA)
  analdati[[roii_pheno]] = y_outliers_removed0$cleaned_data[[roii_pheno]]
  analdati[['y']] = inormal(analdati[[roii_pheno]])# already scaled
  
  covariate_names = covariate_names[which(!psych::describe(subset(analdati,select=covariate_names))$sd==0)]
  mod0.roii = paste0(c('age.c','age.c2',setdiff(covariate_names,c('age','sex'))),
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

get_adj_BMI = function(d,covariate_names){
  analdati = subset(d,select=c('IID','BMI',covariate_names)) %>%
    mutate(age.c = scale(age)[,1]/2) %>% mutate(age.c2=age.c^2) %>%
    mutate(logBMI = log(BMI))#for removing outliers
  
  # remove outliers based on bmi
  x_outliers_removed0 = remove.outliers_grubbs(analdati,varnames="logBMI")
  print(x_outliers_removed0$counts.NA)
  analdati[["logBMI"]] = x_outliers_removed0$cleaned_data[["logBMI"]]
  analdati[['x']] = inormal(analdati[["logBMI"]])
  covariate_names = covariate_names[which(!psych::describe(subset(analdati,select=covariate_names))$sd==0)]
  
  mod0 = paste0(c('age.c','age.c2',setdiff(covariate_names,c('age','sex',"ICV"))),collapse = " + ")
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


# ======= TO-DO: remove all the model-covariates except for sex =======

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

# starting -------------------------------------------------------------------- 
cov_additional= c('current_smoking','hypertension','type2diabetes')
roi_ctx_measurements = c(paste(rois,"volume",sep="_"),
                         paste(rois,"area",sep="_"),
                         paste(rois,"thickness",sep="_"))
#********************* NEED TO CHANGE BEFORE DISTRIBUTION *********************#
#
# filter_i0 = 12; filter_i1=17; filter_is = c(filter_i0,filter_i1)
filter_is = c(12,13)
filter_is = which(df_filter$group_name == 'age <65 years & BMI >=35kg/m2')
filter_is = c(filter_is,which(df_filter$group_name == 'age >=65 years & BMI >=35kg/m2'))
#
#******************************************************************************#

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
    sink(messages, type="message")
    sink(messages, type="output")
    
    adj_ctx = do.call(cbind,lapply(roi_ctx_measurements,get_adj_ctx,
                                   d=d,
                                   covariate_names = cov_ctx))
    adj_BMI = get_adj_BMI(d=d, covariate_names = cov_BMI)
    
    analdat = data.frame( subset(d,select = c("IID",std_covariate_names)),adj_ctx,adj_BMI )
    analdat = get_subset("analdat",filter_expression = filter_expression_i)
    
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
        # APOE4-specific BMI  effect
        smooth_plot_sex_combined[[roii]] <- 
          ggplot(analdati,aes(x=BMI.adj,y=y,color=E4_status,fill=E4_status)) + 
          geom_smooth(alpha=0.2,method = 'gam', formula = y~s(x,bs='cs',k=4)) + 
          # geom_point(alpha=0.1) + 
          ggtitle(unlist(str_split(roii,"_"))[2]) +
          scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
          scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) + xlab("")
        
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
        assoc_res_adj_roii.E4.carrier$Nmax = nrow(analdati %>% filter(E4_status == "E4-carrier"))
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
        assoc_res_adj_roii.E4.noncarrier$Nmax = nrow(analdati %>% filter(E4_status == "E4-noncarrier"))
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
        assoc_res_adj_roii$N = nrow(analdati)
        assoc_res_adj_roii$E4_status = 'APOE4-by-BMI interaction'
        
        # sex-combined BMI-effect
        res.roii = rbind(assoc_res_adj_roii.E4.carrier,
                         assoc_res_adj_roii.E4.noncarrier,
                         assoc_res_adj_roii)
        res.roii$roi = roii
        assoc_res_adj = rbind(assoc_res_adj,res.roii)
      }
      
      ## plotting smooth curves roi-specific BMI-ctx_phenos
      sm_plot_dir = file.path(outdir,paste("smooth_plots_adj",model,"_",filter_group,sep=''))
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
        analdati = subset(analdat,select=c("IID",roii,'logBMI','ICV.scaled','E4_status','age.c','age.c2','sex',
                                           setdiff(c(cov_ctx,cov_BMI),c("age","sex","BMI","ICV")))) 
        
        analdati[['y']] <- inormal(analdati[[paste(roii,sep=".")]])
        analdati[['BMI']] <- inormal(analdati[['logBMI']])
        analdati[['ICV']] <- inormal(analdati[['ICV.scaled']])
        covariate_names = covariate_names[which(!psych::describe(subset(analdati,select=covariate_names))$sd==0)]
        
        mod0 = paste0(c('age.c','age.c2',covariate_names),collapse = " + ")
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
        assoc_res_adj_roii$Nmax = nrow(analdati)
        assoc_res_adj_roii$roi = roii
        assoc_res_adj_roii$E4_status = 'APOE4-by-BMI interaction'
        assoc_res_adj_roii$mod = mod0
        
        rownames(assoc_res_adj_roii)=NULL
        assoc_res_cov = rbind(assoc_res_cov,assoc_res_adj_roii)
      }
      
      # formatting results ----------------------------------------------------------
      res1 = subset(assoc_res_adj,term %in% c("BMI.adj","E4_status:BMI.adj"))
      names(res1)[2:4] = c("Estimate","SE","P")
      res1$method = 'm1'
      res1$term = str_remove(res1$term,"[.]adj")
      
      res2 = assoc_res_cov
      names(res2)[2:4] = c("Estimate","SE","P")
      res2$method = 'm2'
      
      # creating outputs ------------------------------------------------------------
      # table containing both sets of the results
      roi_assoc_res_two_methods = merge(subset(res1,E4_status=="APOE4-by-BMI interaction"),
                                        res2, by=c('roi',"term"))
      names(roi_assoc_res_two_methods) = str_replace(names(roi_assoc_res_two_methods),"[.]x","_reg_out_covs")
      names(roi_assoc_res_two_methods) = str_replace(names(roi_assoc_res_two_methods),"[.]y","_adj_covs_in_mod")
      
      # cleaning-up -----------------------------------------------------------------
      save(
        res1,res2,roi_assoc_res_two_methods,model,group_name,filter_expression_i,filter_group,
        file=file.path(outdir,paste("roi_assoc_res_two_methods_",model,"_",filter_group,".Rdata",sep=''))
      )
      rm(res1,res2,res_two_methods,adj_ctxTH,adj_BMI,analdat,analdati,group_name,filter_expression_i)
    }#if(min(table(analdat$E4_status)) >(30+length(cov_ctx)))
    
    #ending message ---------------------------------------------------------------
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
