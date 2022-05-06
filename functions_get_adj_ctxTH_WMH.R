get_adj_ctxTH = function(roii,d,covariate_names){
  d = d %>% mutate(age.c = scale(age)[,1],
                   age.c2 = (scale(age)/2)^2)# for stability (Gelman)
  
  cov_required = c('age','sex','hypertension','current_smoking','type2diabetes')
  cov_study_specific = setdiff(names(covariate_names),cov_required)
  
  analdati = subset(d,select=c('WMH',names(covariate_names),'age.c', 'age.c2',roii)) 
  analdati[['y']] = rntransform(analdati[[roii]])
  
  analdati = analdati %>% mutate(y = ifelse(abs(scale(y)[,1]) >4, NA, y))
  mod0.roii = paste0(c('age.c','age.c2','hypertension','current_smoking','type2diabetes',cov_study_specific),collapse = " + ")
  mod0.roii = paste("y ~ sex*(",mod0.roii,")",sep="")
  
  # age, sex and other covariates: 
  fit0.roii = lm(as.formula(mod0.roii),data=analdati,na.action=na.exclude)
  analdati[[paste(roii,"adj",sep=".")]] = resid(fit0.roii)
  head(analdati)
  
  ret = subset(analdati,select=paste(roii,"adj",sep="."))
  ret
}

get_adj_WMH = function(d,covariate_names){
  d = d %>% mutate(age.c = scale(age)[,1],
                   age.c2 = (scale(age)/2)^2)# for stability (Gelman)
  
  cov_required = c('age','sex','hypertension','current_smoking','type2diabetes',"ICVorBrainVolume")
  cov_study_specific = setdiff(names(covariate_names),cov_required)
  
  analdati = subset(d,select=c('WMH',names(covariate_names),'age.c', 'age.c2')) %>%
    mutate(ICVorBrainVolume = scale(ICVorBrainVolume)[,1]/2)
  
  analdati[['x']] = rntransform(analdati[["WMH"]])
  
  analdati = analdati %>% mutate(y = ifelse(abs(scale(x)[,1]) >4, NA, x))
  mod0 = paste0(c('age.c','age.c2','hypertension','current_smoking','type2diabetes','ICVorBrainVolume',cov_study_specific),collapse = " + ")
  mod0 = paste("x ~ sex*(",mod0,")",sep="")
  
  # age, sex and other covariates: 
  fit0 = lm(as.formula(mod0),data=analdati,na.action=na.exclude)
  analdati[[paste("WMH","adj",sep=".")]] = resid(fit0)
  head(analdati)
  
  ret = subset(analdati,select=paste("WMH","adj",sep="."))
  ret
}
