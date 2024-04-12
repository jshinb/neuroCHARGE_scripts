#*****************************************************************************#
#
# Step 2: Descriptive statistics: table and plots
#
#*****************************************************************************#
op <- options(warn=1)
logfile=file.path(outdir,"2_Descriptive_statistics.log"); messages=file(logfile, open="wt")
sink(messages, type="message")
sink(messages, type="output")

#1. Descriptive statistics ------------------------------------------------------
cat("\n2. Obtaining descriptive statistics \n") 

desc_All = psych::describe(subset(d,select=-IID) %>% dplyr::select(where(is.numeric)),IQR = T)
desc_FemalesMales = psych::describeBy(subset(d,select=-IID)%>% dplyr::select(where(is.numeric)),group=d$sex,IQR = T)
desc_APOE4status = psych::describeBy(subset(d,select=-IID)%>% dplyr::select(where(is.numeric)),group=d$E4_status,IQR = T)
desc_AgeGroup = psych::describeBy(subset(d,select=-IID)%>% dplyr::select(where(is.numeric)),group=d$age_group,IQR = T)
desc_BMIGroup = psych::describeBy(subset(d,select=-IID)%>% dplyr::select(where(is.numeric)),group=d$BMI_group,IQR = T)
save(desc_All,
     desc_FemalesMales,
     desc_APOE4status,
     desc_AgeGroup,
     desc_BMIGroup,
     file=file.path(outdir,"descriptive_tables.RData"))

# categorical variables:
descriptive_stats_discrete <- function(x,data){
  ret = table(data[[x]])
  ret = data.frame(ret,prop.table(ret))
  ret = subset(ret,select=-Var1.1)
  names(ret) = c('Value',"Count","Proportion")
  ret
}#descriptive_stats_discrete

cat_vars = setdiff(names(d),names(d %>% dplyr::select_if(is.numeric)))

desc_All_discrete = lapply(cat_vars,descriptive_stats_discrete,data=subset(d,select = cat_vars))
names(desc_All_discrete) = cat_vars

desc_Female_discrete = lapply(cat_vars,descriptive_stats_discrete,data=subset(d,select = cat_vars, sex=="F"))
names(desc_Female_discrete) = cat_vars

desc_Male_discrete = lapply(cat_vars,descriptive_stats_discrete,data=subset(d,select = cat_vars, sex=="M"))
names(desc_Male_discrete) = cat_vars

desc_E4carrier_discrete = lapply(cat_vars,descriptive_stats_discrete,data=subset(d,select = cat_vars, E4_status=="E4-carrier"))
names(desc_E4carrier_discrete) = cat_vars

desc_E4noncarrier_discrete = lapply(cat_vars,descriptive_stats_discrete,data=subset(d,select = cat_vars, E4_status=="E4-noncarrier"))
names(desc_E4noncarrier_discrete) = cat_vars

desc_age_gt_65_discrete = lapply(cat_vars,descriptive_stats_discrete,
                                data=subset(d,select = cat_vars, age_group=="age_gt_65"))#age>=65
names(desc_age_gt_65_discrete) = cat_vars

desc_age_lt_65_discrete = lapply(cat_vars,descriptive_stats_discrete,
                                data=subset(d,select = cat_vars, age_group=="age_lt_65"))#age<65
names(desc_age_lt_65_discrete) = cat_vars
##
desc_BMI_gt_25_discrete = lapply(cat_vars,descriptive_stats_discrete,
                                data=subset(d,select = cat_vars, BMI_group=="BMI_gt_25"))#BMI>=25
names(desc_BMI_gt_25_discrete) = cat_vars

desc_BMI_lt_25_discrete = lapply(cat_vars,descriptive_stats_discrete,
                                data=subset(d,select = cat_vars, BMI_group=="BMI_lt_25"))#BMI<25
names(desc_BMI_lt_25_discrete) = cat_vars
##
save(desc_All_discrete,desc_Female_discrete,desc_Male_discrete,
     desc_E4carrier_discrete,desc_E4noncarrier_discrete,
     desc_age_gt_65_discrete,desc_age_lt_65_discrete,
     desc_BMI_gt_25_discrete,desc_BMI_lt_25_discrete,
     file=file.path(outdir,"descriptive_discrete_tables.RData"))

#2. correlation matrices and their plots for 'all' variables --------------------
corm = cor(subset(d,select=setdiff(names(d),c("IID",cat_vars))),use="p")
cormF = cor(subset(d,sex=="F",select=setdiff(names(d),c("IID",cat_vars))),use="p")
cormM = cor(subset(d,sex=="M",select=setdiff(names(d),c("IID",cat_vars))),use="p")

# sex-combined
p = ggcorrplot::ggcorrplot(corm,hc.order = T,hc.method = 'ward.D2',tl.cex = 8)
png(file.path(outdir,"correlation_pairwise_brainphenos_plots.png"),
    width=16, height=16, units="in", res=300)
print(p)
dev.off()
save(corm, file=file.path(outdir,"correlation_matrix_brain_pheno_All.RData"))

# female
p = ggcorrplot::ggcorrplot(cormF,hc.order = T,hc.method = 'ward.D2',tl.cex = 8)
png(file.path(outdir,"correlation_pairwise_brainphenos_plots_Female.png"),
    width=16, height=16, units="in", res=300)
print(p)
dev.off()
save(cormF, file=file.path(outdir,"correlation_matrix_brain_pheno_Female.RData"))

# male
p = ggcorrplot::ggcorrplot(cormM,hc.order = T,hc.method = 'ward.D2',tl.cex = 8)
png(file.path(outdir,"correlation_pairwise_brainphenos_plots_Male.png"),
    width=16, height=16, units="in", res=300)
print(p)
dev.off()
save(cormM, file=file.path(outdir,"correlation_matrix_brain_pheno_Male.RData"))

#3. correlation matrices and their plots for 'all' variables --------------------
corm_E4carrier = cor(subset(d,E4_status=="E4-carrier",select=setdiff(names(d),c("IID",cat_vars))),use="p")
corm_E4noncarrier = cor(subset(d,E4_status == "E4-noncarrier",select=setdiff(names(d),c("IID",cat_vars))),use="p")

p = ggcorrplot::ggcorrplot(corm_E4carrier,hc.order = T,hc.method = 'ward.D2',tl.cex = 8)
png(file.path(outdir,"correlation_pairwise_brainphenos_plots_E4carrier.png"),
    width=16, height=16, units="in", res=300)
print(p)
dev.off()
save(corm_E4carrier, file=file.path(outdir,"correlation_matrix_brain_pheno_E4carrier.RData"))

p = ggcorrplot::ggcorrplot(corm_E4noncarrier,hc.order = T,hc.method = 'ward.D2',tl.cex = 8)
png(file.path(outdir,"correlation_pairwise_brainphenos_plots_E4noncarrier.png"),
    width=16, height=16, units="in", res=300)
print(p)
dev.off()
save(corm_E4noncarrier, file=file.path(outdir,"correlation_matrix_brain_pheno_E4noncarrier.RData"))

n_age_lt_65_E4carrier = nrow(d %>% dplyr::filter(age <65,E4_status=="E4-carrier"))
n_age_gt_65_E4carrier = nrow(d %>% dplyr::filter(age >=65,E4_status=="E4-carrier"))
if(n_age_gt_65_E4carrier-length(std_covariate_names) >=30){
  corm_age_gt_65 = cor(subset(d,age_group=="age_gt_65",select=setdiff(names(d),c("IID",cat_vars))),use="p")
  p = ggcorrplot::ggcorrplot(corm_age_gt_65,hc.order = T,hc.method = 'ward.D2',tl.cex = 8)
  png(file.path(outdir,"correlation_pairwise_brainphenos_plots_age_gt_65.png"),
      width=16, height=16, units="in", res=300)
  print(p)
  dev.off()
  save(corm_age_gt_65, file=file.path(outdir,"correlation_matrix_brain_pheno_age_gt_65.RData"))
}
run_age_gt_65 = n_age_gt_65_E4carrier-length(std_covariate_names) >=30
run_age_lt_65 = n_age_lt_65_E4carrier-length(std_covariate_names) >=30
if(n_age_lt_65_E4carrier-length(std_covariate_names) >=30){
    corm_age_lt_65 = cor(subset(d,age_group=="age_lt_65",select=setdiff(names(d),c("IID",cat_vars))),use="p")
    p = ggcorrplot::ggcorrplot(corm_age_lt_65,hc.order = T,hc.method = 'ward.D2',tl.cex = 8)
    png(file.path(outdir,"correlation_pairwise_brainphenos_plots_age_lt_65.png"),
        width=16, height=16, units="in", res=300)
    print(p)
    dev.off()
    save(corm_age_lt_65, file=file.path(outdir,"correlation_matrix_brain_pheno_age_lt_65.RData"))
  }
# run age-stratified analyses only if there are at least (30+number of covariates) 
# in all the E4-carrier groups
run_age_stratified_analyses = n_age_lt_65_E4carrier-length(std_covariate_names) >=30 & n_age_gt_65_E4carrier-length(std_covariate_names) >=30

n_BMI_lt_25_E4carrier = nrow(d %>% dplyr::filter(BMI <25,E4_status=="E4-carrier"))
n_BMI_gt_25_E4carrier = nrow(d %>% dplyr::filter(BMI >=25,E4_status=="E4-carrier"))
run_BMI_gt_25 = n_BMI_gt_25_E4carrier-length(std_covariate_names) >=30
run_BMI_lt_25 = n_BMI_lt_25_E4carrier-length(std_covariate_names) >=30

if(run_BMI_gt_25){
  corm_BMI_gt_25 = cor(subset(d,BMI_group=="BMI_gt_25",select=setdiff(names(d),c("IID",cat_vars))),use="p")
  p = ggcorrplot::ggcorrplot(corm_BMI_gt_25,hc.order = T,hc.method = 'ward.D2',tl.cex = 8)
  png(file.path(outdir,"correlation_pairwise_brainphenos_plots_BMI_gt_25.png"),
      width=16, height=16, units="in", res=300)
  print(p)
  dev.off()
  save(corm_BMI_gt_25, file=file.path(outdir,"correlation_matrix_brain_pheno_BMI_gt_25.RData"))
}
if(run_BMI_lt_25){
  corm_BMI_lt_25 = cor(subset(d,BMI_group=="BMI_lt_25",select=setdiff(names(d),c("IID",cat_vars))),use="p")
  p = ggcorrplot::ggcorrplot(corm_BMI_lt_25,hc.order = T,hc.method = 'ward.D2',tl.cex = 8)
  png(file.path(outdir,"correlation_pairwise_brainphenos_plots_BMI_lt_25.png"),
      width=16, height=16, units="in", res=300)
  print(p)
  dev.off()
  save(corm_BMI_lt_25, file=file.path(outdir,"correlation_matrix_brain_pheno_BMI_lt_25.RData"))
}
# run BMI-stratified analyses only if there are at least (30+number of covariates) 
# in all the E4-carrier groups
run_BMI_stratified_analyses = run_BMI_gt_25 & run_BMI_lt_25

# sex-combined
p = ggcorrplot::ggcorrplot(corm,hc.order = T,hc.method = 'ward.D2',tl.cex = 8)
png(file.path(outdir,"correlation_pairwise_brainphenos_plots.png"),
    width=16, height=16, units="in", res=300)
print(p)
dev.off()
save(corm, file=file.path(outdir,"correlation_matrix_brain_pheno_All.RData"))
# female
p = ggcorrplot::ggcorrplot(cormF,hc.order = T,hc.method = 'ward.D2',tl.cex = 8)
png(file.path(outdir,"correlation_pairwise_brainphenos_plots_Female.png"),
    width=16, height=16, units="in", res=300)
print(p)
dev.off()
save(cormF, file=file.path(outdir,"correlation_matrix_brain_pheno_Female.RData"))
# male
p = ggcorrplot::ggcorrplot(cormM,hc.order = T,hc.method = 'ward.D2',tl.cex = 8)
png(file.path(outdir,"correlation_pairwise_brainphenos_plots_Male.png"),
    width=16, height=16, units="in", res=300)
print(p)
dev.off()
save(cormM, file=file.path(outdir,"correlation_matrix_brain_pheno_Male.RData"))

#4. scatter plots with correlation ----------------------------------------------
dA=subset(d,select=std_covariate_names)
png(file.path(outdir,"correlation_pairwise_scatter_All.png"),
    width=11,height=11,units="in",res=300)
pairs.panels(data.frame(dA,stringsAsFactors = F), 
             method = "pearson", # correlation method
             hist.col = "#BEBEBE33",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
dev.off()

png(file.path(outdir,"correlation_pairwise_scatter_Female.png"),
    width=11,height=11,units="in",res=300)
pairs.panels(data.frame(subset(dA,sex=="F"),stringsAsFactors = F), 
             method = "pearson", # correlation method
             hist.col = "#8B000033",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
dev.off()

png(file.path(outdir,"correlation_pairwise_scatter_Male.png"),
    width=11,height=11,units="in",res=300)
pairs.panels(data.frame(subset(dA,sex=="M"),stringsAsFactors = F), 
             method = "pearson", # correlation method
             hist.col = "#00008B33",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
dev.off()


#5. characteristic tables -------------------------------------------------------
# 
tab_vars = setdiff(std_covariate_names,"IID")
tableOne_All <- CreateTableOne(vars = tab_vars,
                               #strata = "trt", 
                               data = d, 
                               factorVars = cat_vars)
save(tableOne_All,file=file.path(outdir,"characteristic_table.RData"))

tableOne_sex_str <- CreateTableOne(vars = tab_vars,
                                   strata = "sex", 
                                   data = d, 
                                   factorVars = cat_vars)
save(tableOne_sex_str,file=file.path(outdir,"characteristic_table_sex_stratified.RData"))

tableOne_APOE4_str <- CreateTableOne(vars = tab_vars,
                                   strata = "E4_status", 
                                   data = d, 
                                   factorVars = cat_vars)
save(tableOne_APOE4_str,file=file.path(outdir,"characteristic_table_APOE4_stratified.RData"))

if(n_age_gt_65_E4carrier>0 | n_age_lt_65_E4carrier>0){
  tableOne_age_str <- CreateTableOne(vars = tab_vars,
                                     strata = "age_group", 
                                     data = d, 
                                     factorVars = cat_vars)
  save(tableOne_age_str,file=file.path(outdir,"characteristic_table_age_stratified.RData"))
}

if(n_BMI_gt_25_E4carrier>0 | n_BMI_lt_25_E4carrier>0){
  tableOne_BMI_str <- CreateTableOne(vars = tab_vars,
                                     strata = "BMI_group", 
                                     data = d, 
                                     factorVars = cat_vars)
  save(tableOne_BMI_str,file=file.path(outdir,"characteristic_table_BMI_stratified.RData"))
}

#6. family size ---------------------------------------------------------------
if(any(ls()=="FID")){
  if(!is.na(FID)){
    fam.size.list = vector(mode='list',length=3) 
    
    fam.size.list[[1]] = table(table(d[["FID"]]))
    
    fam.size_E4_status = table(d[["FID"]],d$E4_status)
    fam.size.list[[2]] = table(paste(fam.size_E4_status[,"E4-noncarrier"],fam.size_E4_status[,"E4-carrier"],sep="_"))
    
    fam.size_BMI = table(d[["FID"]],d$BMI_group)
    fam.size.list[[3]] = table(paste(fam.size_BMI[,"BMI_lt_25"],fam.size_BMI[,"BMI_gt_25"],sep="_"))
    
    names(fam.size.list) = c('all',#'nE4_noncarrier_vs_carrier',
                             paste0(colnames(table(d[["FID"]],d$E4_status)),collapse="_vs_"),
                             paste0(colnames(table(d[["FID"]],d$BMI_group)),collapse="_vs_"))
    save(fam.size.list,file=file.path(outdir,"fam_size_tables.RData"))
  }
}

sink()
closeAllConnections()
print(readLines(logfile))

cat("\nFinishing obtaining descriptive statistics.\n") 
options(op)
