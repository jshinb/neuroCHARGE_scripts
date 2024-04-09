#*****************************************************************************#
#
# 3a: Under age<65 years:
# Examine associations between BMI vs. ctx-phenotypes across the 34 regions
# adjusting for age, sex, ICV (for brain phenotypes) and any 
# cohort-specific covariates in the participants <65 years old.
#
model='baseModel'
age_group = "age_lt_65"
#*****************************************************************************#
op <- options(warn=1)
cat("\n3a. Start obtaining association estimates for BMI vs. roi-ctx phenotypes for the base model (age <65 years).\n")
logfile=file.path(
  outdir,
  paste("3a_",age_group,"_roiAssociation_statistics_",model,".log",sep='')
  ); messages=file(logfile, open="wt")

sink(messages, type="message")
sink(messages, type="output")

# starting -------------------------------------------------------------------- 
rois = c("bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", 
         "cuneus", "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal", 
         "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal", 
         "lingual", "medialorbitofrontal", "middletemporal", "parahippocampal", 
         "paracentral", "parsopercularis", "parsorbitalis", "parstriangularis", 
         "pericalcarine", "postcentral", "posteriorcingulate", "precentral", 
         "precuneus", "rostralanteriorcingulate", "rostralmiddlefrontal", 
         "superiorfrontal", "superiorparietal", "superiortemporal", "supramarginal", 
         "frontalpole", "temporalpole", "transversetemporal", "insula")  

roi_ctx_measurements = c(paste(rois,"volume",sep="_"),paste(rois,"area",sep="_"),paste(rois,"thickness",sep="_"))
cov_additional= c('current_smoking','hypertension','type2diabetes')
cov_ctx = setdiff(std_covariate_names,c("BMI","age_BMI","E4_status",cov_additional))
cov_BMI = setdiff(std_covariate_names,c("BMI","age_BMI","E4_status","ICV",cov_additional))

# approach1 -------------------------------------------------------------------
adj_ctx = do.call(cbind,lapply(roi_ctx_measurements,get_adj_ctx,
                               d=d %>% filter(age<65),
                                 covariate_names = cov_ctx))
adj_BMI = get_adj_BMI(d=d %>% filter(age<65), covariate_names = cov_BMI)

## create empty lists for outputs 
smooth_plot_sex_combined <- smooth_plot_sex_stratified <- vector(mode="list",length=length(roi_ctx_measurements))
names(smooth_plot_sex_combined) <- names(smooth_plot_sex_stratified) <- roi_ctx_measurements
assoc_res_adj = c()

## run fitting base models
for(roii in roi_ctx_measurements){
  analdati = data.frame(subset(d %>% filter(age<65),
                               select = c("IID",std_covariate_names)),
                        adj_BMI,
                        subset(adj_ctx,select=paste(roii,"adj",sep=".")))
  
  analdati[['y']] <- analdati[[paste(roii,"adj",sep=".")]]
  analdati = analdati %>% filter(!is.na(E4_status))
  # APOE4-specific BMI  effect
  smooth_plot_sex_combined[[roii]] <- analdati %>% 
    ggplot(aes(x=BMI.adj,y=y,color=E4_status,fill=E4_status)) + 
    geom_smooth(alpha=0.2,method = 'gam') + 
    # geom_point(alpha=0.1) + 
    ggtitle(roii) + 
    scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
    scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))
    
  smooth_plot_sex_stratified[[roii]] <- analdati %>% 
    ggplot(aes(x=BMI.adj,y=y,color=E4_status,fill=E4_status,linetype=sex)) + 
    geom_smooth(alpha=0.2,method = 'gam') + 
    # geom_point(alpha=0.1) + 
    ggtitle(roii) + 
    scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) +
    scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))
  
  ## E4-carrier
  fit.assoc.E4.carrier = lm(y~BMI.adj+sex,data=analdati,subset=E4_status=="E4-carrier", 
                            na.action = na.exclude)
  excl_col = which(colnames(summary(fit.assoc.E4.carrier)$coef)=="t value")
  assoc_res_adj_roii.E4.carrier = summary(fit.assoc.E4.carrier)$coef["BMI.adj",-excl_col,drop=F]
  terms = rownames(assoc_res_adj_roii.E4.carrier)
  rownames(assoc_res_adj_roii.E4.carrier)=NULL
  assoc_res_adj_roii.E4.carrier = data.frame(term=terms,assoc_res_adj_roii.E4.carrier)
  assoc_res_adj_roii.E4.carrier$N = nobs(fit.assoc.E4.carrier)
  assoc_res_adj_roii.E4.carrier$E4_status ="E4-carrier"
  
  ## E4-noncarrier
  fit.assoc.E4.noncarrier = lm(y~BMI.adj+sex,data=analdati,subset=E4_status=="E4-noncarrier", 
                               na.action = na.exclude)
  excl_col = which(colnames(summary(fit.assoc.E4.noncarrier)$coef)=="t value")
  assoc_res_adj_roii.E4.noncarrier = summary(fit.assoc.E4.noncarrier)$coef["BMI.adj",-excl_col,drop=F]
  terms = rownames(assoc_res_adj_roii.E4.noncarrier)
  rownames(assoc_res_adj_roii.E4.noncarrier)=NULL
  assoc_res_adj_roii.E4.noncarrier = data.frame(term=terms,assoc_res_adj_roii.E4.noncarrier)
  assoc_res_adj_roii.E4.noncarrier$N = nobs(fit.assoc.E4.noncarrier)
  assoc_res_adj_roii.E4.noncarrier$E4_status ="E4-noncarrier"
  
  # E4_status-by-BMI interaction effect
  fit.assoc = lm(y~E4_status*BMI.adj+sex,data=analdati,na.action = na.exclude)
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

## plotting smooth curves roi-specific BMI-ctx_phenos
if(file.exists(file.path(outdir,"smooth_plots_adjBase_age_lt_65"))){
  unlink(file.path(outdir,"smooth_plots_adjBase_age_lt_65"),recursive = T)
}

dir.create(file.path(outdir,"smooth_plots_adjBase_age_lt_65"))
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
    ylab("adjusted for base covariates") 
  
  sm_roi_common_ylim = sm_roi & 
    coord_cartesian(ylim=range(p_ranges_y))
  
  pfile = file.path(outdir,'smooth_plots_adjBase_age_lt_65',paste("sm_plot_adjBase_",roi,"%03d_age_lt_65.png",sep=""))
  png(file=pfile,width=11,height=8,units='in',res=300)
  print(sm_roi)
  print(sm_roi_common_ylim)
  dev.off()
}

# approach2 -------------------------------------------------------------------
analdat = d %>% mutate(age.c = scale(age)[,1],
                 age.c2 = (scale(age)/2)^2) %>% # for stability (Gelman)
  mutate(logBMI = log(BMI), ICV.scaled = scale(ICV)[,1]/2)
## removing outliers of BMI and ICV
d_cleaned = remove.outliers_grubbs(analdat,varnames = c('logBMI',"ICV.scaled"))
analdat = d_cleaned$cleaned_data
print(d_cleaned$counts.NA)
rm(d_cleaned)

assoc_res_cov = c()
for(roii in roi_ctx_measurements){
  analdati = subset(analdat%>%filter(age<65),select=c("IID",roii,'logBMI','ICV.scaled','E4_status','age.c','age.c2','sex',
                               setdiff(c(cov_ctx,cov_BMI),c("age","sex","BMI","ICV")))) 

  analdati[['y']] <- inormal(analdati[[paste(roii,sep=".")]])
  analdati[['BMI']] <- inormal(analdati[['logBMI']])
  analdati[['ICV']] <- inormal(analdati[['ICV.scaled']])
  
  mod0 = paste0(c('age.c','age.c2',setdiff(c(cov_ctx,cov_BMI),c('age','sex'))),collapse = " + ")
  mod0 = paste("y ~ E4_status*BMI + sex*(",mod0,")",sep="")
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

adj_ctx_base_age_lt_65 = adj_ctx
adj_BMI_base_age_lt_65 = adj_BMI

# cleaning-up -----------------------------------------------------------------
save(
  res1,res2,roi_assoc_res_two_methods,model,age_group,
  file=file.path(outdir,paste("roi_assoc_res_two_methods_",model,"_",age_group,".Rdata",sep=''))
)

rm(p,res1,res2,res_two_methods,adj_ctxTH,adj_BMI,analdat,analdati)

sink()
closeAllConnections()
print(readLines(logfile))

#ending message ---------------------------------------------------------------
cat("\nFinished obtaining association estimates for BMI vs. roi-ctx phenotypes under base model.\n")
options(op)
