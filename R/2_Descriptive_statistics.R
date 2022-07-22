#*****************************************************************************#
#
# Step 2: Descriptive statistics: table and plots
#
#*****************************************************************************#
# Descriptive statistics ------------------------------------------------------
cat("\n2. Obtaining descriptive statistics \n") 

desc_All = psych::describe(subset(d,select=-IID),IQR = T)
desc_FemalesMales = psych::describeBy(subset(d,select=-IID),group="sex",IQR = T)
save(desc_All,desc_FemalesMales,file=file.path(outdir,"descriptive_tables.RData"))

# categorical variables:
descriptive_stats_discrete <- function(x,data){
  ret = table(data[[x]])
  ret = data.frame(ret,prop.table(ret))
  ret = subset(ret,select=-Var1.1)
  names(ret) = c('Value',"Count","Proportion")
  ret
}#descriptive_stats_discrete

cat_vars = rownames(desc_All)[str_detect(rownames(desc_All),"[*]")]
cat_vars = str_remove(cat_vars,"[*]")
desc_All_discrete = lapply(cat_vars,descriptive_stats_discrete,data=subset(d,select = cat_vars))
names(desc_All_discrete) = cat_vars

desc_Female_discrete = lapply(cat_vars,descriptive_stats_discrete,data=subset(d,select = cat_vars, sex=="F"))
names(desc_Female_discrete) = cat_vars

desc_Male_discrete = lapply(cat_vars,descriptive_stats_discrete,data=subset(d,select = cat_vars, sex=="M"))
names(desc_Male_discrete) = cat_vars
save(desc_All_discrete,desc_Female_discrete,desc_Male_discrete,
     file=file.path(outdir,"descriptive_discrete_tables.RData"))

# correlation matrices and their plots for 'all' variables --------------------
corm = cor(subset(d,select=setdiff(names(d),c("IID",cat_vars))),use="p")
cormF = cor(subset(d,sex=="F",select=setdiff(names(d),c("IID",cat_vars))),use="p")
cormM = cor(subset(d,sex=="M",select=setdiff(names(d),c("IID",cat_vars))),use="p")

png(file.path(outdir,"correlation_pairwise_brainphenos_plots.png"),
    width=12, height=3.5, units="in", res=300)
par(mfrow=c(1,3))
# sex-combined
corrplot(corm,tl.cex=0.7,tl.col = "darkgrey",
         bg = "White", #tl.srt = 25, 
         title = "\n\n sex-combined \n",
         #addCoef.col = "black", 
         type = "full")# save correlation matrix, too
save(corm, file=file.path(outdir,"correlation_matrix_brain_pheno_All.RData"))
# female
corrplot(cormF,tl.cex=0.7,tl.col = "darkgrey",
         bg = "White", #tl.srt = 25, 
         title = "\n\n female \n",
         #addCoef.col = "black", 
         type = "full")# save correlation matrix, too
save(cormF, file=file.path(outdir,"correlation_matrix_brain_pheno_Female.RData"))
# male
corrplot(cormM,tl.cex=0.7,tl.col = "darkgrey",
         bg = "White", #tl.srt = 25, 
         title = "\n\n male \n",
         #addCoef.col = "black", 
         type = "full")# save correlation matrix, too
save(cormM, file=file.path(outdir,"correlation_matrix_brain_pheno_Male.RData"))
dev.off()


# scatter plots with correlation ----------------------------------------------
dA=subset(d,select=setdiff(c(names(covariate_names),"logWMH","ICVorBrainVolume"),c('IID')))
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


# characteristic tables -------------------------------------------------------
# 
tab_vars = setdiff(c(names(sel_colnames),'logWMH'),"IID")
tableOne_All <- CreateTableOne(vars = tab_vars,
                               #strata = "trt", 
                               data = d, 
                               factorVars = cat_vars)

tableOne_sex_str <- CreateTableOne(vars = tab_vars,
                                   strata = "sex", 
                                   data = d, 
                                   factorVars = cat_vars)

save(tableOne_All,file=file.path(outdir,"characteristic_table.RData"))
save(tableOne_sex_str,file=file.path(outdir,"characteristic_table_stratified.RData"))

cat("\nFinishing obtaining descriptive statistics.\n") 
cat("\n# -------------------------------------------------------------------------------------- #\n",
    file=file.path(outdir,input_specification_file), append=T)
cat("Warnings:\n",file=file.path(outdir,input_specification_file), append=T)
capture.output(summary(warnings()),
               file=file.path(outdir,input_specification_file), append=T)
