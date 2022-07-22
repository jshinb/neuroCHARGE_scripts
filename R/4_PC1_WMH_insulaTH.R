#*****************************************************************************#
#
# 4: PC1 of cov-adjusted WMH and insulaTH 
#
#*****************************************************************************#
  cat("\n4. Starting calculation of PC1 for base model covariate-adjusted WMH and insular thickness.\n") 

pca_analdat = data.frame(adj_WMH_base,subset(adj_ctxTH_base,select='insula.adj'),
                         row.names = d$IID)
# out1 - pca plot
myPCA.na.omit = PCA(na.omit(pca_analdat),graph=F)
p = plot.PCA(myPCA.na.omit,choix = 'var')
p = p + theme(panel.grid.major = element_blank(),
          plot.title=element_text(size=14, color="darkblue"),
          axis.title = element_text(size=10, color=gray(0.4),hjust=0.5))

pdf(file.path(outdir, 'PCA_graph_WMH_InsuralTH.pdf'), width=6.5, height=5)
print(p)
dev.off()

GWAS.pheno = merge(data.frame(IID=rownames(pca_analdat),pca_analdat),
                   data.frame(IID=rownames(myPCA.na.omit$ind$coord),PC1=myPCA.na.omit$ind$coord[,1]),
                   by="IID",all.x=T)
GWAS.pheno = GWAS.pheno %>% mutate(IID = as.character(IID))

# creating histograms of each variable: insular TH, WMH, and PC1
plot_d = merge(subset(d,select=c(IID,sex)) %>% mutate(IID = as.character(IID)),
               GWAS.pheno,
               by="IID")
hist.WMH = plot_d %>% 
  ggplot(aes(x=WMH.adj,fill=sex)) + 
  geom_histogram( color="#e9ecef", alpha=0.4, position = 'identity' ) + 
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="")+
  ggtitle("covariate-adjusted WMH") + 
  theme(legend.position = 'none')

hist.ctx= plot_d %>% 
  ggplot(aes(x=insula.adj,fill=sex)) + 
  geom_histogram( color="#e9ecef", alpha=0.4, position = 'identity' ) + 
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="") + 
  ggtitle("covariate-adjusted insular thickness") + 
  theme(legend.position = 'none')

hist.PC1= plot_d %>% 
  ggplot(aes(x=PC1,fill=sex)) + 
  geom_histogram( color="#e9ecef", alpha=0.4, position = 'identity' ) + 
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="") + 
  ggtitle("PC1")
p_hist = hist.WMH + hist.ctx + hist.PC1

png(file.path(outdir,"hist_adjusted_WMH_insularTH_PC1.png"),
    width=8.5,height=3,units='in',res=300)
print(p_hist)
dev.off()

# pca - eigenvalues, pct of variance explained, and loadings
PCA_eigenV =  data.frame(myPCA.na.omit$eig)
PCA_loadings = data.frame(myPCA.na.omit$var$coord)
rownames(PCA_loadings) = paste("loading",rownames(PCA_loadings),sep=".")
PCA_res = data.frame(PCA_eigenV,t(PCA_loadings))
save(PCA_res,file=file.path(outdir,"PCA_res.Rdata"))

# sign: 
if(sign(PCA_loadings["loading.insula.adj","Dim.1"])==1){
  cat("\n ** Flipping sign of PC1 to be correlated -vely with insularTH.**\n",
      file=file.path(outdir,input_specification_file),append=T)
  plot_d$PC1 = (-1)*plot_d$PC1
}
corm_PC = cor(subset(plot_d,select=c(WMH.adj,insula.adj,PC1)),use="p")
save(corm_PC,file=file.path(outdir,"corm_PC.Rdata"))


# save the derived phenotypes: Use these derived phenotypes for GWAS analyses 
# with study-specific covariates
gwas_pheno_dir = paste(cohort_name,ancestry,"GWAS_file",sep="_")
dir.create(gwas_pheno_dir)
write_tsv(GWAS.pheno,file.path(gwas_pheno_dir,"covariate_adjusted_phenotypes_for_GWAS.tsv"))
capture.output(describe(GWAS.pheno),
               file=file.path(outdir,input_specification_file), append=T)

cat("\nFinishing calculation of PC1 for covariate-adjusted WMH and insular thickness.\n") 

cat("\n# -------------------------------------------------------------------------------------- #\n",
    file=file.path(outdir,input_specification_file), append=T)
cat("Warnings:\n",file=file.path(outdir,input_specification_file), append=T)
capture.output(summary(warnings()),
               file=file.path(outdir,input_specification_file), append=T)
