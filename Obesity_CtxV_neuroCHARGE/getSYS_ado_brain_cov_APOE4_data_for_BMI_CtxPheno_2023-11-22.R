pacman::p_load(tidyverse, data.table, psych, 
               patchwork, GGally,corrplot,hrbrthemes,
               tableone,FactoMineR,factoextra,
               ppcor, lsmeans, multcomp, ModelMetrics,
               #caret, 
               gridExtra, Hmisc, pastecs, testit)
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
wd = '/Users/jshin/Library/CloudStorage/OneDrive-SickKids/neuroCHARGE_WMH_ctxTH_age/WMH_ctxTH/Example_SYS'
setwd(wd)
sys_data_dir ='/Users/jshin/Library/CloudStorage/OneDrive-SickKids/SYS_Database'
