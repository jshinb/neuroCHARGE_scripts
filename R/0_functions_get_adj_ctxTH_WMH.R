#' Inverse normal transformation (INT)
#'
#' https://www.biostars.org/p/80597/
#' See the supplement of Yang et al. Nature 2012. 
#'
#' @example
#' x1 <- 5:1
#' inormal(x1)
#'
#' x2 <- c(NA, 5:1, 5:1, NA)
#' inormal(x2)
inormal <- function(x)
{
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

get_adj_ctxTH = function(roii,d,covariate_names){
  d = d %>% mutate(age.c = scale(age)[,1],
                   age.c2 = (scale(age)/2)^2)# for stability (Gelman)
  
  cov_required = c('age','sex','hypertension','current_smoking','type2diabetes')
  cov_study_specific = setdiff(names(covariate_names),cov_required)
  
  analdati = subset(d,select=c('WMH',names(covariate_names),'age.c', 'age.c2',roii)) 
  analdati[['y']] = inormal(analdati[[roii]])
  
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
  
  analdati[['x']] = inormal(analdati[["WMH"]])
  
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
