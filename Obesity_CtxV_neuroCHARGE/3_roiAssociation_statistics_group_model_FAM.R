#*****************************************************************************#
#
# 3 FAM:
# Examine associations between BMI vs. ctx-phenotype across the 34 regions
# adjusting for age, sex, ICV (for brain phenotype) and any 
# cohort-specific covariates in the participants according to the filter grouping.
#
models=c('base','full')
filter_groups = c("all","age_gt_65","age_lt_65","BMI_gt_25","BMI_lt_25")
filter_group_names = c("all","age >=65 years","age <65 years", "BMI >=25kg/m2", "BMI <25kg/m2")
#
#*****************************************************************************#
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

for(model in models){
  if(model=="base"){
    cov_BMI = setdiff(std_covariate_names,c("BMI","age_BMI","E4_status","ICV",cov_additional))
    cov_ctx = setdiff(std_covariate_names,c("BMI","age_BMI","E4_status",cov_additional))
  }else{#full
    cov_ctx = setdiff(std_covariate_names,c("BMI","age_BMI","E4_status"))
    cov_BMI = setdiff(std_covariate_names,c("BMI","age_BMI","E4_status","ICV"))
  }
  for(filter_group in filter_groups){
    gi = which(filter_groups == filter_group)
    msg = paste("\n3. Start estimating associations between BMI vs. roi-ctx phenotype under ",
                model," model in the participants with ",filter_group_names[gi],".\n",sep='')
    cat(msg)
    logfile=file.path(outdir,paste("3_",filter_group,"_roiAssociation_statistics_",model,".log",sep='')); messages=file(logfile, open="wt")
    sink(messages, type="message")
    sink(messages, type="output")
    
    if(filter_group=="age_gt_65"){
      analdat = d %>% filter(age >=65)
    }else if(filter_group=="age_lt_65"){
      analdat = d %>% filter(age <65)
    }else if(filter_group=="BMI_gt_25"){
      analdat = d %>% filter(BMI >=25)
    }else if(filter_group=="BMI_lt_25"){
      analdat = d %>% filter(BMI <25)
    }else{#'all'
      analdat = d
    }
    
    if(min(table(analdat$E4_status)) >(30+length(cov_ctx))){
      # approach1 -------------------------------------------------------------------
      ## create empty lists for outputs 
      smooth_plot_sex_combined <- smooth_plot_sex_stratified <- vector(mode="list",length=length(roi_ctx_measurements))
      names(smooth_plot_sex_combined) <- names(smooth_plot_sex_stratified) <- roi_ctx_measurements
      assoc_res_adj = c()
      
      adj_ctx = do.call(cbind,lapply(roi_ctx_measurements,get_adj_ctx,
                                     d=analdat,
                                     covariate_names = cov_ctx))
      adj_BMI = get_adj_BMI(d=analdat, covariate_names = cov_BMI)
      
      
      ## run fitting base models
      glob_ylab <- paste('adjusted brain phenotype')
      glob_xlab <- paste("adjusted BMI")
      for(roii in roi_ctx_measurements){
        analdati = data.frame(subset(analdat,select = c("IID","FID",std_covariate_names)),
                              adj_BMI,subset(adj_ctx,select=paste(roii,"adj",sep=".")))
        
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
        
        ## E4-carrier
        fit.assoc.E4.carrier = lmer(y~BMI.adj+sex+(1|FID),data=analdati,subset=E4_status=="E4-carrier", 
                                  na.action = na.exclude)
        excl_col = which(colnames(summary(fit.assoc.E4.carrier)$coef)=="t value")
        assoc_res_adj_roii.E4.carrier = summary(fit.assoc.E4.carrier)$coef["BMI.adj",-excl_col,drop=F]
        terms = rownames(assoc_res_adj_roii.E4.carrier)
        rownames(assoc_res_adj_roii.E4.carrier)=NULL
        assoc_res_adj_roii.E4.carrier = data.frame(term=terms,assoc_res_adj_roii.E4.carrier)
        assoc_res_adj_roii.E4.carrier$N = nobs(fit.assoc.E4.carrier)
        assoc_res_adj_roii.E4.carrier$Nfam = summary(fit.assoc.E4.carrier)$ngrps
        assoc_res_adj_roii.E4.carrier$SD_fam= attr(summary(fit.assoc.E4.carrier)$varcor$FID,"stddev")
        assoc_res_adj_roii.E4.carrier$E4_status ="E4-carrier"
        
        ## E4-noncarrier
        fit.assoc.E4.noncarrier = lmer(y~BMI.adj+sex+(1|FID),data=analdati,subset=E4_status=="E4-noncarrier", 
                                     na.action = na.exclude)
        excl_col = which(colnames(summary(fit.assoc.E4.noncarrier)$coef)=="t value")
        assoc_res_adj_roii.E4.noncarrier = summary(fit.assoc.E4.noncarrier)$coef["BMI.adj",-excl_col,drop=F]
        terms = rownames(assoc_res_adj_roii.E4.noncarrier)
        rownames(assoc_res_adj_roii.E4.noncarrier)=NULL
        assoc_res_adj_roii.E4.noncarrier = data.frame(term=terms,assoc_res_adj_roii.E4.noncarrier)
        assoc_res_adj_roii.E4.noncarrier$N = nobs(fit.assoc.E4.noncarrier)
        assoc_res_adj_roii.E4.noncarrier$Nfam = summary(fit.assoc.E4.noncarrier)$ngrps
        assoc_res_adj_roii.E4.noncarrier$SD_fam= attr(summary(fit.assoc.E4.noncarrier)$varcor$FID,"stddev")
        assoc_res_adj_roii.E4.noncarrier$E4_status ="E4-noncarrier"
        
        # E4_status-by-BMI interaction effect
        fit.assoc = lmer(y~E4_status*BMI.adj+sex + (1|FID),data=analdati,na.action = na.exclude)
        excl_col = which(colnames(summary(fit.assoc)$coef)=="t value")
        assoc_res_adj_roii = summary(fit.assoc)$coef[,-excl_col,drop=F]
        terms = rownames(assoc_res_adj_roii)
        rownames(assoc_res_adj_roii)=NULL
        assoc_res_adj_roii = data.frame(term=terms,assoc_res_adj_roii)
        assoc_res_adj_roii$N = nobs(fit.assoc)
        assoc_res_adj_roii$Nfam = summary(fit.assoc)$ngrps
        assoc_res_adj_roii$SD_fam= attr(summary(fit.assoc)$varcor$FID,"stddev")
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
      
      # approach2 -------------------------------------------------------------------
      analdat = d %>% mutate(age.c = scale(age)[,1],
                             age.c2 = (scale(age)/2)^2) %>% # for stability (Gelman)
        mutate(logBMI = log(BMI), ICV.scaled = scale(ICV)[,1]/2)
      if(filter_group=="age_gt_65"){
        analdat = analdat %>% filter(age >=65)
      }else if(filter_group=="age_lt_65"){
        analdat = analdat %>% filter(age <65)
      }else if(filter_group=="BMI_gt_25"){
        analdat = analdat %>% filter(BMI >=25)
      }else if(filter_group=="BMI_lt_25"){
        analdat = analdat %>% filter(BMI <25)
      }else{# 'all'
        analdat = analdat
      }
      ## removing outliers of BMI and ICV
      d_cleaned = remove.outliers_grubbs(analdat,varnames = c('logBMI',"ICV.scaled"))
      analdat = d_cleaned$cleaned_data
      print(d_cleaned$counts.NA)
      rm(d_cleaned)
      
      assoc_res_cov = c()
      mods = c()
      for(roii in roi_ctx_measurements){
        analdati = subset(analdat,select=c("IID","FID",roii,'logBMI','ICV.scaled','E4_status','age.c','age.c2','sex',
                                           setdiff(c(cov_ctx,cov_BMI),c("age","sex","BMI","ICV")))) 
        
        analdati[['y']] <- inormal(analdati[[paste(roii,sep=".")]])
        analdati[['BMI']] <- inormal(analdati[['logBMI']])
        analdati[['ICV']] <- inormal(analdati[['ICV.scaled']])
        
        mod0 = paste0(c('age.c','age.c2',setdiff(c(cov_ctx,cov_BMI),c('age','sex'))),collapse = " + ")
        mod0 = paste("y ~ E4_status*BMI + sex*(",mod0,") + (1|FID)",sep="")
        mods = c(mods,mod0)
        capture.output(mod0,file=file.path(outdir,input_specification_file), append=T)
        
        fit.assoc = lmer(mod0,data=analdati,na.action = na.exclude)
        assoc_res_adj_roii = summary(fit.assoc)$coef
        terms = rownames(assoc_res_adj_roii)
        excl_col = which(colnames(assoc_res_adj_roii)=="t value")
        assoc_res_adj_roii = data.frame(term=terms,assoc_res_adj_roii[,-excl_col])
        assoc_res_adj_roii = subset(assoc_res_adj_roii, term %in% c("BMI","E4_statusE4-carrier:BMI"))
        assoc_res_adj_roii$N = nobs(fit.assoc)
        assoc_res_adj_roii$Nfam = summary(fit.assoc)$ngrps
        assoc_res_adj_roii$SD_fam= attr(summary(fit.assoc)$varcor$FID,"stddev")
        assoc_res_adj_roii$roi = roii
        assoc_res_adj_roii$E4_status = 'APOE4-by-BMI interaction'
        assoc_res_adj_roii$mod = mod0
        
        rownames(assoc_res_adj_roii)=NULL
        assoc_res_cov = rbind(assoc_res_cov,assoc_res_adj_roii)
      }#for(roii in roi_ctx_measurements)
      
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
        res1,res2,roi_assoc_res_two_methods,model,filter_group,
        file=file.path(outdir,paste("roi_assoc_res_two_methods_",model,"_",filter_group,".Rdata",sep=''))
      )
      
      # adj_ctx_base_BMI_lt_25 = adj_ctx
      # adj_BMI_base_BMI_lt_25 = adj_BMI
      rm(p,res1,res2,res_two_methods,adj_ctxTH,adj_BMI,analdat,analdati,filter_group)
    }#if(min(table(analdat$E4_status)) >(30+length(cov_ctx)))
    
    #ending message ---------------------------------------------------------------
    end.msg = paste("\nFinished obtaining association estimates for BMI vs. roi-ctx phenotypes under ",model," model (",filter_group_names[gi],").\n",sep='')
    cat(end.msg)
    
    sink()
    closeAllConnections()
    print(readLines(logfile))
  } #for(filter_group in filter_groups)
  rm(model)
} #for(model in models)
options(op)
