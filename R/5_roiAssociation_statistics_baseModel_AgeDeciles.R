#*****************************************************************************#
#
# 5. examine associations between WWH vs. ctxTH across the 34 regions
# adjusting for the basic covariates and cohort-specific covariates.
#
#*****************************************************************************#
cat("\n5. Start obtaining association estimates for WMH vs. roi-ctxTH for the base model.\n")
# starting -------------------------------------------------------------------- 

roi.34 = setdiff(names(brainpheno_names),c("WMH"))
cov_additional= c('current_smoking','hypertension','type2diabetes','bmi')
cov_ctxTH = names(covariate_names)[!names(covariate_names) %in% c('ICVorBrainVolume',cov_additional)]
cov_WMH = names(covariate_names)[!names(covariate_names) %in% cov_additional]

d = subset(d,!is.na(age)) %>% mutate(age.c = scale(age)[,1],
                 age.c2 = (scale(age)/2)^2)# for stability (Gelman)
b2 = quantile(d$age,probs = seq(from=0,to=1,length.out=11))
d = d %>% 
  mutate(age_cat2 = cut(age,breaks=b2,include.lowest=T))

# -----------------------------------------------------------------------------

cov_required = covariate_names[unique(c(cov_ctxTH,cov_WMH))]
assoc_res_cov = c()
age_cats = unique(d$age_cat2)
for(age_cati in age_cats){
  di = subset(d,age_cat2 == age_cati);
  for(roii in roi.34){
    analdati = subset(di,select=c(roii,'WMH','age.c','age.c2','sex',
                                 setdiff(names(cov_required),c("age","sex"))))
    analdati[['y']] <- inormal(analdati[[paste(roii,sep=".")]])
    analdati[['WMH']] <- inormal(analdati[["WMH"]])
    analdati <- analdati %>% 
      mutate(ICVorBrainVolume = scale(ICVorBrainVolume)[,1]/2)
    #sex-combined-interaction
    mod0 = paste0(c('WMH','age.c','age.c2',setdiff(names(cov_required),c('age','sex'))),collapse = " + ")
    mod0 = paste("y ~ sex*(",mod0,")",sep="")
    capture.output(mod0,file=file.path(outdir,input_specification_file), append=T)
    
    fit.assoc = lm(mod0,data=analdati,na.action = na.exclude)
    assoc_res_adj_roii = summary(fit.assoc)$coef
    terms = rownames(assoc_res_adj_roii)
    assoc_res_adj_roii = data.frame(term=terms,assoc_res_adj_roii[,-3])
    assoc_res_adj_roii = subset(assoc_res_adj_roii, term %in% c("WMH","sexM:WMH"))
    assoc_res_adj_roii$N = nobs(fit.assoc)
    assoc_res_adj_roii$roi = roii
    assoc_res_adj_roii$sex = 'sex-by-WMH interaction'
    assoc_res_adj_roii$mod = mod0
    assoc_res_adj_roii$age_cat = age_cati    
    assoc_res_cov = rbind(assoc_res_cov,assoc_res_adj_roii)
    
    # sex-stratified
    mod0 = paste0(c('WMH','age.c','age.c2',setdiff(names(cov_required),c('age','sex'))),collapse = " + ")
    mod0 = paste("y ~ ",mod0,sep="")
    capture.output(mod0,file=file.path(outdir,input_specification_file), append=T)
    
    ##female
    fit.assoc = lm(mod0,data=analdati,subset=sex=="F",na.action = na.exclude)
    assoc_res_adj_roii = summary(fit.assoc)$coef
    terms = rownames(assoc_res_adj_roii)
    assoc_res_adj_roii = data.frame(term=terms,assoc_res_adj_roii[,-3])
    assoc_res_adj_roii = subset(assoc_res_adj_roii, term =="WMH")
    assoc_res_adj_roii$N = nobs(fit.assoc)
    assoc_res_adj_roii$roi = roii
    assoc_res_adj_roii$sex = 'female'
    assoc_res_adj_roii$mod = mod0
    assoc_res_adj_roii$age_cat = age_cati    
    assoc_res_cov = rbind(assoc_res_cov,assoc_res_adj_roii)
    
    ##male
    fit.assoc = lm(mod0,data=analdati,subset=sex=="M",na.action = na.exclude)
    assoc_res_adj_roii = summary(fit.assoc)$coef
    terms = rownames(assoc_res_adj_roii)
    assoc_res_adj_roii = data.frame(term=terms,assoc_res_adj_roii[,-3])
    assoc_res_adj_roii = subset(assoc_res_adj_roii, term =="WMH")
    assoc_res_adj_roii$N = nobs(fit.assoc)
    assoc_res_adj_roii$roi = roii
    assoc_res_adj_roii$sex = 'male'
    assoc_res_adj_roii$mod = mod0
    assoc_res_adj_roii$age_cat = age_cati    
    assoc_res_cov = rbind(assoc_res_cov,assoc_res_adj_roii)

    # sex-combined
    mod0 = paste0(c('WMH','age.c','age.c2',setdiff(names(cov_required),c('age'))),collapse = " + ")
    mod0 = paste("y ~ ",mod0,sep="")
    capture.output(mod0,file=file.path(outdir,input_specification_file), append=T)
    
    fit.assoc = lm(mod0,data=analdati,na.action = na.exclude)
    assoc_res_adj_roii = summary(fit.assoc)$coef
    terms = rownames(assoc_res_adj_roii)
    assoc_res_adj_roii = data.frame(term=terms,assoc_res_adj_roii[,-3])
    assoc_res_adj_roii = subset(assoc_res_adj_roii, term =="WMH")
    assoc_res_adj_roii$N = nobs(fit.assoc)
    assoc_res_adj_roii$roi = roii
    assoc_res_adj_roii$sex = 'sex-combined'
    assoc_res_adj_roii$mod = mod0
    assoc_res_adj_roii$age_cat = age_cati    
    
    assoc_res_cov = rbind(assoc_res_cov,assoc_res_adj_roii)
  }
}

# formatting results ----------------------------------------------------------
res2 = assoc_res_cov
names(res2)[2:4] = c("Estimate","SE","P")
res2$method = 'm2'
write_tsv(res2,file.path(outdir,"roi_assoc_res_baseModel_age_deciles.tsv"))

# creating outputs ------------------------------------------------------------
age_cat_levels = rev(names(table(d$age_cat2)))

L=length(age_cat_levels)
COLS <- colorRampPalette(c("black","orange"))(L)
names(COLS) <- age_cat_levels

res2_all = res2 %>% 
  mutate(age_cat = factor(age_cat,levels = age_cat_levels)) %>%
  arrange(age_cat)
res2 = subset(res2_all,term=="WMH" & sex!="sex-by-WMH interaction") %>%
  mutate(sex = factor(sex,levels = (c('sex-combined','female','male'))))

# roi.34.ordered = subset(res2,sex=="sex-combined" & age_cat==age_cat_levels[1]) %>% 
#   arrange(Estimate)
# roi.34.ordered = roi.34.ordered$roi
roi.34.ordered=c("insula", "superiortemporal", "fusiform", "parsopercularis", 
                 "lateralorbitofrontal", "inferiortemporal", "medialorbitofrontal", 
                 "rostralanteriorcingulate", "supramarginal", "middletemporal", 
                 "parsorbitalis", "MCT", "parahippocampal", "caudalanteriorcingulate", 
                 "rostralmiddlefrontal", "isthmuscingulate", "lingual", "transversetemporal", 
                 "temporalpole", "entorhinal", "posteriorcingulate", "bankssts", 
                 "parstriangularis", "lateraloccipital", "frontalpole", "precentral", 
                 "superiorfrontal", "inferiorparietal", "pericalcarine", "paracentral", 
                 "caudalmiddlefrontal", "cuneus", "precuneus", "postcentral", 
                 "superiorparietal")#from base model
res2 = res2 %>% 
  mutate(roi=factor(roi,levels=rev(roi.34.ordered))) %>% 
  mutate(U95CI = Estimate + 1.96*SE,L95CI = Estimate - 1.96*SE) %>%
  arrange(roi)

pos <- position_dodge(width=0)#(width=0.75)#
p <- res2 %>% 
  ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,color=age_cat,group=age_cat)) + 
  facet_grid(cols = vars(sex)) + 
  geom_point(size=1, stroke=0.5, position=pos,alpha=0.8) +
  geom_errorbar(position=pos, width=0.25,size=0.5,alpha=0.5) + 
  geom_hline(yintercept=0, size=0.4) +
  geom_path(aes(group=age_cat),position=pos,alpha=0.5) +
  geom_hline(yintercept = 0)+
  xlab(NULL) + 
  coord_flip() +
  scale_color_manual(values = COLS) +
  ylim((range(res2$L95CI,res2$U95CI)))+#only when coord_flip
  ggtitle("WMH vs. cortical thickness associations") + 
  guides(col = guide_legend(reverse=T))

pdf(file.path(outdir,"forest_plot_assoc_WMH_roi_ctxTH_baseModel_deciles.pdf"),
    width=18,height=6.5)
print(p + theme_bw())
dev.off()

#ending message ---------------------------------------------------------------
cat("\nFinished obtaining association estimates for WMH vs. roi-ctxTH under base model in age-based deciles.\n")

cat("\n# -------------------------------------------------------------------------------------- #\n",
    file=file.path(outdir,input_specification_file), append=T)
cat("Warnings:\n",file=file.path(outdir,input_specification_file), append=T)
capture.output(summary(warnings()),
               file=file.path(outdir,input_specification_file), append=T)

rm(p,res2,res2_all)
