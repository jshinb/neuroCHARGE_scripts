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
rois = c("bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", 
         "cuneus", "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal", 
         "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal", 
         "lingual", "medialorbitofrontal", "middletemporal", "parahippocampal", 
         "paracentral", "parsopercularis", "parsorbitalis", "parstriangularis", 
         "pericalcarine", "postcentral", "posteriorcingulate", "precentral", 
         "precuneus", "rostralanteriorcingulate", "rostralmiddlefrontal", 
         "superiorfrontal", "superiorparietal", "superiortemporal", "supramarginal", 
         "frontalpole", "temporalpole", "transversetemporal", "insula")  

remove.outliers_grubbs <- function(dat, varnames){
  require(outliers)
  count.na0 <- count.na1 <- c()
  
  for(varname in varnames){
    x = dat[[varname]]
    count.na0 = c(count.na0,sum(is.na(x)))
    keep.going <- T
    while(keep.going){
      test <- grubbs.test(x,opposite = T)
      print(test)
      if(test$p.value<0.05){
        if(str_detect(test$alternative,"lowest")){
          cat("-------------------------------------------------------------------")
          cat(paste('\nReplacing the lowest value of \'',varname,'\' with NA.\n',sep=''))
          cat("-------------------------------------------------------------------\n")
          x[which.min(x)] <- NA
        }else if(str_detect(test$alternative,"highest")){
          cat("-------------------------------------------------------------------")
          cat(paste('\nReplacing the highest value of \'',varname,'\' with NA.\n',sep=''))
          cat("-------------------------------------------------------------------")
          x[which.max(x)] <- NA
        }
      }
      test2 <- grubbs.test(x,opposite = F)
      print(test2)
      if(test2$p.value<0.05){
        if(str_detect(test2$alternative,"lowest")){
          cat("-------------------------------------------------------------------")
          cat(paste('\nReplacing the lowest value of \'',varname,'\' with NA.\n',sep=''))
          cat("-------------------------------------------------------------------")
          x[which.min(x)] <- NA
        }else if(str_detect(test2$alternative,"highest")){
          cat("-------------------------------------------------------------------")
          cat(paste('\nReplacing the highest value of \'',varname,'\' with NA.\n',sep=''))
          cat("-------------------------------------------------------------------")
          x[which.max(x)] <- NA
        }
      }
      
      keep.going=test$p.value<0.05 | test2$p.value<0.05
    }
    count.na1 = c(count.na1,sum(is.na(x)))
    dat[[varname]] <- x
  }
  counts.NA = data.frame(var=varnames,before=count.na0,after=count.na1)
  ret = list(counts.NA = counts.NA, cleaned_data=dat)
  ret
}

get_VolSums <- function(data,ID_col){#data -> cortical volume data 
  rois = c("bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", 
           "cuneus", "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal", 
           "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal", 
           "lingual", "medialorbitofrontal", "middletemporal", "parahippocampal", 
           "paracentral", "parsopercularis", "parsorbitalis", "parstriangularis", 
           "pericalcarine", "postcentral", "posteriorcingulate", "precentral", 
           "precuneus", "rostralanteriorcingulate", "rostralmiddlefrontal", 
           "superiorfrontal", "superiorparietal", "superiortemporal", "supramarginal", 
           "frontalpole", "temporalpole", "transversetemporal", "insula")  
  SumVol = subset(data,select=ID_col)
  for(roi in rois){
    Lh = paste("lh",roi,"volume",sep="_")
    Rh = paste("rh",roi,"volume",sep="_")
    SumVol_roi = (data[[Lh]] + data[[Rh]])
    SumVol = cbind(SumVol,SumVol_roi) #34 rois
  }
  SumVol = data.frame(SumVol) %>% mutate(TotalVol = apply(data.frame(SumVol)[,-1],1,sum,na.rm=T))
  names(SumVol) = c(ID_col,paste(rois,"volume",sep="_"),"total_volume")
  SumVol
}

get_ThicknessMeans <- function(data,ID_col){#data -> cortical volume data 
  rois =c("bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", 
          "cuneus", "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal", 
          "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal", 
          "lingual", "medialorbitofrontal", "middletemporal", "parahippocampal", 
          "paracentral", "parsopercularis", "parsorbitalis", "parstriangularis", 
          "pericalcarine", "postcentral", "posteriorcingulate", "precentral", 
          "precuneus", "rostralanteriorcingulate", "rostralmiddlefrontal", 
          "superiorfrontal", "superiorparietal", "superiortemporal", "supramarginal", 
          "frontalpole", "temporalpole", "transversetemporal", "insula")
  MeanTH = subset(data,select=ID_col)
  for(roi in rois){
    Lh = paste("lh",roi,"thickness",sep="_")
    Rh = paste("rh",roi,"thickness",sep="_")
    MeanTH_roi = (data[[Lh]] + data[[Rh]])/2
    MeanTH = cbind(MeanTH,MeanTH_roi) #34 rois
  }
  Global = (data[['lh_MeanThickness_thickness']] + 
              data[['rh_MeanThickness_thickness']])/2
  MeanTH = data.frame(MeanTH, Global)
  names(MeanTH) = c(ID_col,paste(rois,"thickness",sep="_"),"global_mean_thickness")
  MeanTH
}

get_AreaSums <- function(data,ID_col){#data -> cortical volume data
  rois =c("bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal",
          "cuneus", "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal",
          "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal",
          "lingual", "medialorbitofrontal", "middletemporal", "parahippocampal",
          "paracentral", "parsopercularis", "parsorbitalis", "parstriangularis",
          "pericalcarine", "postcentral", "posteriorcingulate", "precentral",
          "precuneus", "rostralanteriorcingulate", "rostralmiddlefrontal",
          "superiorfrontal", "superiorparietal", "superiortemporal", "supramarginal",
          "frontalpole", "temporalpole", "transversetemporal", "insula")
  
  Sum = subset(data,select=ID_col)
  for(roi in rois){
    Lh = paste("lh",roi,"area",sep="_")
    Rh = paste("rh",roi,"area",sep="_")
    Sum_roi = (data[[Lh]] + data[[Rh]])
    Sum = cbind(Sum,Sum_roi) #34 rois
  }
  
  Total = data[['rh_WhiteSurfArea_area']] + data[['lh_WhiteSurfArea_area']]
  Sum = data.frame(Sum, Total)
  names(Sum) = c(ID_col,paste(rois,"area",sep="_"),"total_area")
  Sum
}

inormal <- function(x)
{
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

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

get_subset = function(data_name, filter_expression){
  if(!is.na(filter_expression)){
    execute_string = paste("analdat <- ",data_name,"%>% filter(",filter_expression,")",sep='')
    eval(parse(text=execute_string))
  }else{
    execute_string = paste("analdat <- ",data_name,sep='')
    eval(parse(text=execute_string))
  }#else
  analdat
}

models=c('base','full')
