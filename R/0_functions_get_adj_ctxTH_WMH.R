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
get_means <- function(data,ID_col){#data -> cortical thickness data 
  rois = fread('freesurfer_34rois.txt',header=F)$V1
  avgTH = subset(data,select=ID_col)
  for(roi in rois){
    Lh = paste("lh",roi,"thickness",sep="_")
    Rh = paste("rh",roi,"thickness",sep="_")
    avgTH_roi = (data[[Lh]] + data[[Rh]])/2
    avgTH = cbind(avgTH,avgTH_roi) #34 rois
  }
  L_SA = 'lh_WhiteSurfArea_area'
  R_SA = 'rh_WhiteSurfArea_area'
  L_TH = 'lh_MeanThickness_thickness'
  R_TH = 'rh_MeanThickness_thickness'
  totalSA = data[[L_SA]] + data[[R_SA]]
  global_avgTH = data[[L_TH]]*(data[[L_SA]]/totalSA) + data[[R_TH]]*(data[[R_SA]]/totalSA)
  avgTH = data.frame(avgTH, global_avgTH)
  names(avgTH) = c(ID_col,rois,"MCT")
  avgTH
}

inormal <- function(x)
{
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

get_adj_ctxTH = function(roii,d,covariate_names){
  
  d = d %>% mutate(age.c = scale(age)[,1],
                   age.c2 = (scale(age)/2)^2)# for stability (Gelman)
  
  analdati = subset(d,select=c('WMH',names(covariate_names),'age.c', 'age.c2',roii)) 
  analdati[['y']] = inormal(analdati[[roii]])
  
  analdati = analdati %>% mutate(y = ifelse(abs(scale(y)[,1]) >4, NA, y))
  mod0.roii = paste0(c('age.c','age.c2',setdiff(names(covariate_names),c('age','sex'))),collapse = " + ")
  mod0.roii = paste("y ~ sex*(",mod0.roii,")",sep="")
  
  # age, sex and other covariates: 
  fit0.roii = lm(as.formula(mod0.roii),data=analdati,na.action=na.exclude)
  analdati[[paste(roii,"adj",sep=".")]] = resid(fit0.roii)
  head(analdati)
  paste(mod0.roii)

  ret = subset(analdati,select=paste(roii,"adj",sep="."))
  ret
}

get_adj_WMH = function(d,covariate_names){
  analdati = subset(d,select=c('IID','WMH',names(covariate_names))) %>%
    mutate(age.c = scale(age)[,1],ICVorBrainVolume = scale(ICVorBrainVolume)[,1]/2)

  analdati[['x']] = inormal(analdati[["WMH"]])
  # remove outliers based on WMH
  analdati = analdati %>% mutate(y = ifelse(abs(scale(x)[,1]) >4, NA, x))
  mod0 = paste0(c('age.c',setdiff(names(covariate_names),c('age','sex'))),collapse = " + ")
  mod0 = paste("x ~ sex*(",mod0,")",sep="")
  
  # age, sex and other covariates: 
  fit0 = lm(as.formula(mod0),data=analdati,na.action=na.exclude)
  analdati[[paste("WMH","adj",sep=".")]] = resid(fit0)
  head(analdati)
  paste(mod0)
  
  ret = subset(analdati,select=paste("WMH","adj",sep="."))
  ret
}
