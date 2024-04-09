# ******************** SPECIFY INPUT FILE HERE ******************************** 
#
#
input_specification_file = 'cohort_specific_inputs_UKBB.txt'
#
#
# *****************************************************************************
#-----------------------------------------------------------------------------#

# install/load libraries ------------------------------------------------------
cat("Prep: installing/loading libraries\n")

# install pacman package for 'installing' and 'loading' packages
if(!is.element("pacman",installed.packages()[,1])){
  install.packages("pacman")
}

pacman::p_load(tidyverse, tidyselect, data.table, readxl, psych, outliers,
               patchwork, ggcorrplot, hrbrthemes, gridExtra, ggseg,
               tableone, FactoMineR, factoextra,
               ppcor, lsmeans, multcomp, ModelMetrics,
               # GGally,caret, 
               Hmisc, pastecs, testit,
               futile.logger,tryCatchLog)# for error/warning messages


if(length(find.package("ggsegDefaultExtra", quiet=TRUE, verbose=TRUE))==0){
  remotes::install_github("LCBC-UiO/ggsegDefaultExtra")
}
library(ggsegDefaultExtra)

#source the input-specification file ------------------------------------------
# need to think about what to do
op <- options(nwarnings = 10000, warn=0)
setwd("./") #current directory
source(input_specification_file)

outdir=paste(cohort_name,ancestry,"outputs",sep="_")
# # remove the output directory if a previous folder exists
# if(file.exists(outdir)){
#   unlink(outdir,recursive = T)
# }
dir.create(outdir)
file.copy(input_specification_file,outdir,overwrite = T)

# define functions ------------------------------------------------------------
cat("Prep: sourcing functions\n")
source("0_functions_get_adj_ctxVol_BMI.R")
df_filter = read_xlsx('df_filter.xlsx')

bmi.strata.lelves = c("BMI <20kg/m2",
                      "20<= BMI <25kg/m2", 
                      "25<= BMI <30kg/m2", 
                      "30<= BMI <35kg/m2", 
                      "35<= BMI <40kg/m2", 
                      "BMI >=40kg/m2")

filter_i0 = 12; filter_i1=17
res1_all = c()
for (model_i in c('base','full')){
  for(filter_i in filter_i0:filter_i1){
    group_name = df_filter$group_name[filter_i]
    filter_group = df_filter$group_label[filter_i]
    filter_expression_i = df_filter$expression[filter_i]
    load(file.path(outdir,paste("roi_assoc_res_two_methods_",model_i,"_",filter_group,".Rdata",sep='')))
    res1_all = rbind(res1_all,res1 %>% mutate(model = model_i, group = group_name))
  }
}
print(range(res1_all$Estimate))# -0.5831485  0.7298184

# plotting --------------------------------------------------------------------
plotd = res1_all %>% 
  mutate(U95CI = Estimate + 1.96*SE,L95CI = Estimate - 1.96*SE) %>%
  mutate(pheno = str_split(res1_all$roi,"_",simplify=T)[,2]) %>%
  mutate(roi = str_split(res1_all$roi,"_",simplify=T)[,1]) %>%
  filter(E4_status != 'APOE4-by-BMI interaction') %>% 
  mutate(E4_status = factor(E4_status,levels=c("E4-noncarrier","E4-carrier"))) %>%
  mutate(group = factor(group, levels = (bmi.strata.lelves)))
print(range(plotd$L95CI,plotd$U95CI))#-1.057374  1.270331
yrange_forest = 1.05*c(-1,1)*1.270331#0.4032839

cleanup = theme(panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_line(colour = "grey"),
                axis.line = element_line(colour = "black"),
                panel.border = element_rect(colour = "grey", fill=NA, linewidth=0.5))

## sample sizes ---------------------------------------------------------------
plotd %>% ggplot(aes(x=group, y=N, color=E4_status, shape=model)) + geom_point() + 
  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))

## forest plot ----------------------------------------------------------------
plot.alpha = 0.7
pos <- position_dodge(width=0.75)#

for (model_plot in c("base","full")){
  for (pheno_plot in c("volume","area","thickness") ){
    p <- plotd %>% filter(model==model_plot, pheno==pheno_plot) %>% arrange(roi,group) %>%
      ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,
                 #color=E4_status, alpha=group, 
                 color = group,
                 group=group, shape=group)) + 
      geom_point(size=1, stroke=0.5, position=pos,alpha=plot.alpha) +
      geom_errorbar(position=pos, width=0.25,size=0.5,alpha=plot.alpha) + 
      geom_hline(yintercept=0, size=0.4) +
      facet_grid(E4_status~.)+
      # facet_grid(.~E4_status)+
      # geom_path(aes(group=group),position=pos,alpha=0.5) +
      geom_hline(yintercept = 0, linewidth=0.5)+
      xlab(NULL) + 
      # coord_flip() +
      ylim(yrange_forest) + #only when coord_flip
      # ylim((range(res1_pheno$L95CI,res1_pheno$U95CI)))+#only when coord_flip
      ggtitle(paste("BMI vs. cortical",pheno_plot,"associations")) + 
      # scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) + 
      guides(col = guide_legend(reverse=F,title = element_blank()), 
             shape= guide_legend(reverse=F,title = element_blank()))
    p + theme(legend.position = 'top',axis.text.x = element_text(angle = -30, hjust=0, vjust=0.5)) + cleanup + 
      theme(plot.margin = margin(t=0.15, r=1.5, b=0.15, l=0.15, "cm")) + scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
    
    ggsave(file.path(outdir,paste("forest_plot_by_BMI_classes_",model_plot,"_",pheno_plot,".png",sep='')),
           width=12*0.8,height=8.5*0.8,units='in',dpi=300)
    # width=8.5,height=11,units='in',dpi=300)
  }
}