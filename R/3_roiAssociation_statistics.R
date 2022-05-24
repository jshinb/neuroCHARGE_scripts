#*****************************************************************************#
#
# Step 3: examine associations between WWH vs. ctxTH across the 34 regions
# adjust cortical thickness for age, sex, hypertension, current smoking and any 
# existing covariate-specific covariates 
#
#*****************************************************************************#
cat("\nStep3: Start obtaining association estimates for WMH vs. roi-ctxTH.\n")
# starting -------------------------------------------------------------------- 

roi.34 = setdiff(names(brainpheno_names),c("WMH"))

# approach1 -------------------------------------------------------------------
adj_ctxTH = do.call(cbind,lapply(roi.34,get_adj_ctxTH,d=d,covariate_names = covariate_names))
adj_WMH = get_adj_WMH(d=d,covariate_names = covariate_names)

assoc_res_adj = c()
for(roii in roi.34){
  analdati = data.frame(subset(d,select = c(IID,sex)),
                        adj_WMH,
                        subset(adj_ctxTH,select=paste(roii,"adj",sep=".")))
  analdati[['y']] <- analdati[[paste(roii,"adj",sep=".")]]
  
  # sex-specific WMH  effect
  ## female
  fit.assoc.F = lm(y~WMH.adj,data=analdati,subset=sex=="F", na.action = na.exclude)
  assoc_res_adj_roii.F = summary(fit.assoc.F)$coef["WMH.adj",-3,drop=F]
  terms = rownames(assoc_res_adj_roii.F)
  rownames(assoc_res_adj_roii.F)=NULL
  assoc_res_adj_roii.F = data.frame(term=terms,assoc_res_adj_roii.F)
  assoc_res_adj_roii.F$N = nobs(fit.assoc.F)
  assoc_res_adj_roii.F$sex ="F"
  
  # male
  fit.assoc.M = lm(y~WMH.adj,data=analdati,subset=sex=="M", na.action = na.exclude)
  assoc_res_adj_roii.M = summary(fit.assoc.M)$coef["WMH.adj",-3,drop=F]
  terms = rownames(assoc_res_adj_roii.M)
  rownames(assoc_res_adj_roii.M)=NULL
  assoc_res_adj_roii.M = data.frame(term=terms,assoc_res_adj_roii.M)
  assoc_res_adj_roii.M$N = nobs(fit.assoc.M)
  assoc_res_adj_roii.M$sex ="M"
  
  # sex-combined
  fit.assoc.sex_combined = lm(y~WMH.adj+sex,data=analdati,na.action = na.exclude)
  assoc_res_adj_roii.sex_combined = summary(fit.assoc.sex_combined)$coef["WMH.adj",-3,drop=F]
  terms = rownames(assoc_res_adj_roii.sex_combined)
  rownames(assoc_res_adj_roii.sex_combined)=NULL
  assoc_res_adj_roii.sex_combined = data.frame(term=terms,assoc_res_adj_roii.sex_combined)
  assoc_res_adj_roii.sex_combined$N = nobs(fit.assoc.sex_combined)
  assoc_res_adj_roii.sex_combined$sex ="sex-combined"

  # sex-by-WMH interaction effect
  fit.assoc = lm(y~sex*WMH.adj,data=analdati,na.action = na.exclude)
  assoc_res_adj_roii = summary(fit.assoc)$coef[,-3,drop=F]
  terms = rownames(assoc_res_adj_roii)
  rownames(assoc_res_adj_roii)=NULL
  assoc_res_adj_roii = data.frame(term=terms,assoc_res_adj_roii)
  assoc_res_adj_roii$N = nobs(fit.assoc)
  assoc_res_adj_roii$sex = 'sex-by-WMH interaction'
  
  # sex-combined WMH-effect
  res.roii = rbind(assoc_res_adj_roii.F,
                   assoc_res_adj_roii.M,
                   assoc_res_adj_roii.sex_combined,
                   assoc_res_adj_roii)
  res.roii$roi = roii
  assoc_res_adj = rbind(assoc_res_adj,res.roii)
}

# approach2 -------------------------------------------------------------------
assoc_res_cov = c()
d = d %>% mutate(age.c = scale(age)[,1],
                 age.c2 = (scale(age)/2)^2)# for stability (Gelman)
cov_required = c('age','sex','hypertension','current_smoking','type2diabetes',"ICVorBrainVolume")
cov_study_specific = setdiff(names(covariate_names),cov_required)
for(roii in roi.34){
  analdati = subset(d,select=c(roii,'WMH','ICVorBrainVolume','age.c','age.c2','sex','hypertension','current_smoking','type2diabetes',cov_study_specific))
  analdati[['y']] <- inormal(analdati[[paste(roii,sep=".")]])
  analdati[['WMH']] <- inormal(analdati[["WMH"]])
  analdati <- analdati %>% mutate(ICVorBrainVolume = scale(ICVorBrainVolume)[,1]/2)
  mod0 = paste0(c('WMH','age.c','age.c2','hypertension','current_smoking','type2diabetes','ICVorBrainVolume',cov_study_specific),collapse = " + ")
  mod0 = paste("y ~ sex*(",mod0,")",sep="")
  
  fit.assoc = lm(mod0,data=analdati,na.action = na.exclude)
  assoc_res_adj_roii = summary(fit.assoc)$coef
  terms = rownames(assoc_res_adj_roii)
  assoc_res_adj_roii = data.frame(term=terms,assoc_res_adj_roii[,-3])
  assoc_res_adj_roii = subset(assoc_res_adj_roii, term %in% c("WMH","sexM:WMH"))
  assoc_res_adj_roii$N = nobs(fit.assoc)
  assoc_res_adj_roii$roi = roii
  assoc_res_adj_roii$sex = 'sex-by-WMH interaction'
  assoc_res_cov = rbind(assoc_res_cov,assoc_res_adj_roii)
}

# formatting results ----------------------------------------------------------
res1 = subset(assoc_res_adj,term %in% c("WMH.adj","sexM:WMH.adj"))
names(res1)[2:4] = c("Estimate","SE","P")
res1$method = 'm1'
res1$term = str_remove(res1$term,"[.]adj")

res2 = assoc_res_cov
names(res2)[2:4] = c("Estimate","SE","P")
res2$method = 'm2'


# creating outputs ------------------------------------------------------------
# table containing both sets of the results
res_two_methods = merge(subset(res1,sex=="sex-by-WMH interaction"),
                        res2, by=c('roi',"term"))
roi_assoc_res_two_methods = res_two_methods
save(res1,res2,roi_assoc_res_two_methods,
     file=file.path(outdir,"roi_assoc_res_two_methods.Rdata"))

# plot
beta.p = res_two_methods %>% ggplot(aes(x=Estimate.x,y=Estimate.y,color=term)) +
  geom_point() + geom_abline(intercept = 0, slope = 1) + 
  theme(legend.position = 'none')
SE.p = res_two_methods %>% ggplot(aes(x=SE.x,y=SE.y,color=term)) +
  geom_point() + geom_abline(intercept = 0, slope = 1)+ 
  theme(legend.position = 'none')
Pval.p = res_two_methods %>% ggplot(aes(x=-log10(P.x),y=-log10(P.y),color=term)) +
  geom_point() + geom_abline(intercept = 0, slope = 1)
p_method_comp = beta.p + SE.p + Pval.p
png(file.path(outdir,'roi_assoc_res_two_methods.png'),
    width=8.5,height = 3, units="in", res=300)
print(p_method_comp+ plot_annotation(
  title = 'Comparing two approaches to adjusting for covariates',
  subtitle = 'Method 1 (x-axis): correct first -> examine associations vs. Method2 (y-axis): correct covariates while fitting.'
))
dev.off()

roi.34.ordered = subset(res1,sex=="sex-combined") %>% arrange(Estimate)
roi.34.ordered = roi.34.ordered$roi
res1 = res1 %>% mutate(roi=factor(roi,levels=rev(roi.34.ordered))) %>% 
  mutate(U95CI = Estimate + 1.96*SE,L95CI = Estimate - 1.96*SE) %>%
  arrange(roi)

pos <- position_dodge(width=0.75)#
p <- subset(res1,sex!="sex-by-WMH interaction") %>% 
  ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,color=sex,group=sex)) + 
  geom_point(size=1, stroke=0.5, position=pos,alpha=0.5) +
  geom_errorbar(position=pos, width=0.25,size=0.5,alpha=0.5) + 
  geom_hline(yintercept=0, size=0.4) +
  geom_path(aes(group=sex),alpha=0.5) +
  #scale_color_manual(values=COLS) + 
  geom_hline(yintercept = 0)+
  xlab(NULL) + 
  #ylab("Estimate (adjusted for imaging centre)")+ 
  # theme(legend.position = 'none')+
  # theme(axis.text.x = element_text(angle = 90,hjust =1,vjust=0.5)) + #when not coord_flip() 
  # coord_cartesian(ylim=(range(age_specific_coef$Estimate))) + #when not coord_flip()
  coord_flip() +
  ylim((range(res1$L95CI,res1$U95CI)))+#only when coord_flip
  ggtitle("WMH vs. cortical thickness associations")

pdf(file.path(outdir,"forest_plot_assoc_WMH_roi_ctxTH.pdf"),
    width=6,height=6.5)
print(p)
dev.off()

#ending message ---------------------------------------------------------------
cat("Finished obtaining association estimates for WMH vs. roi-ctxTH.\n")

cat("\n# -------------------------------------------------------------------------------------- #\n",
    file=file.path(outdir,input_specification_file), append=T)
cat("Step3 - Warnings:\n",file=file.path(outdir,input_specification_file), append=T)
capture.output(summary(warnings()),
               file=file.path(outdir,input_specification_file), append=T)

