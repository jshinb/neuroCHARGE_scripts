# ******************** SPECIFY INPUT FILE HERE ******************************** 
#
#
input_specification_file='cohort_specific_inputs_UKBB.txt'
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
               patchwork, ggcorrplot, hrbrthemes, gridExtra, ggseg, viridis,
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
op <- options(nwarnings=10000, warn=0)
setwd("./") #current directory
source(input_specification_file)

outdir=paste(cohort_name,ancestry,"outputs",sep="_")
outdir2=NULL#'adj_sex_only'
# outdir2='adj_sex_only'
# # remove the output directory if a previous folder exists
# if(file.exists(outdir)){
#   unlink(outdir,recursive=T)
# }
dir.create(outdir)
file.copy(input_specification_file,outdir,overwrite=T)

# define functions ------------------------------------------------------------
cat("Prep: sourcing functions\n")
source("0_functions_get_adj_ctxVol_BMI.R")
df_filter=read_xlsx('df_filter.xlsx')
source("add_sample_sizes_2023-12-21.R")

# filter_i0=12; filter_i1=17;filter_is=filter_i0:filter_i1
# bmi.strata.lelves=c("BMI <20kg/m2",
#                       "20<= BMI <25kg/m2", 
#                       "25<= BMI <30kg/m2", 
#                       "30<= BMI <35kg/m2", 
#                       "35<= BMI <40kg/m2", 
#                       "BMI >=40kg/m2")
bmi.strata.lelves=
  # all
  c('BMI <20kg/m2',
  "20<= BMI <25kg/m2",
  "25<= BMI <35kg/m2", 
  "BMI >=35kg/m2",
  # youger
  "age <65 years & BMI <20kg/m2", 
  "age <65 years & 20<= BMI <25kg/m2", 
  "age <65 years & 25<= BMI <35kg/m2", 
  "age <65 years & BMI >=35kg/m2", 
  # older
  "age >=65 years & BMI <20kg/m2", 
  "age >=65 years & 20<= BMI <25kg/m2", 
  "age >=65 years & 25<= BMI <35kg/m2", 
  "age >=65 years & BMI >=35kg/m2")
# sel_bmi_levels=c(1:4); age_group="all"
# sel_bmi_levels=c(5:8); age_group="age_lt_65"
# sel_bmi_levels=c(9:12); age_group="age_gt_65"
# filter_is=(df_filter %>% filter(group_name %in% bmi.strata.lelves[sel_bmi_levels]))$index #c(12,13, 18:21,26,27,32:35)
# a=(df_filter %>% filter(group_name %in% bmi.strata.lelves[sel_bmi_levels]))$group_name
# b=(df_filter %>% filter(group_name %in% bmi.strata.lelves[sel_bmi_levels]))$N_E4_noncarrier
# c=(df_filter %>% filter(group_name %in% bmi.strata.lelves[sel_bmi_levels]))$N_E4_carrier
# write_tsv(data.frame(rbind(a,b,c)), "~/Downloads/tmp.tsv")

filter_is=which(df_filter$group_name == 'age <65 years & BMI >=35kg/m2')
filter_is=c(filter_is,which(df_filter$group_name == 'age >=65 years & BMI >=35kg/m2'))
res1_all=c()
ranges_Estimate=c()
ranges_95CI=c()
# for (model_i in c('base','full')){
for (model_i in c('base')){
  for(filter_i in filter_is){
    sampleN_E4noncarrier=df_filter$N_E4_noncarrier[filter_i]
    sampleN_E4carrier=df_filter$N_E4_carrier[filter_i]
    group_name=df_filter$group_name[filter_i]
    filter_group=df_filter$group_label[filter_i]
    filter_expression_i=df_filter$expression[filter_i]
    if(is.null(outdir2)){
      outf=file.path(outdir,paste("roi_assoc_res_two_methods_",model_i,"_",filter_group,".Rdata",sep=''))
    }else{
      outf=file.path(outdir,outdir2,paste("roi_assoc_res_two_methods_",model_i,"_",filter_group,".Rdata",sep=''))
    }
    if(file.exists(outf)){
      load(outf)
    }else{
      stop(cat("\n",outf," does not exist.\n",sep=''))
    }
    res1_all=rbind(res1_all,res1 %>% mutate(model=model_i, group=group_name, 
                                              Nmax_E4noncarrier=sampleN_E4noncarrier,
                                              Nmax_E4carrier=sampleN_E4carrier,
                                              U95CI=Estimate + 1.96*SE,L95CI=Estimate - 1.96*SE))
    ranges_Estimate=rbind(ranges_Estimate, range(res1_all$Estimate))
    ranges_95CI=rbind(ranges_95CI, range(res1_all$L95CI,res1_all$U95CI))
  }
}

res1_all=res1_all %>% 
  mutate(group=factor(group, levels=bmi.strata.lelves))
data.frame(table(res1_all$group))
view(res1_all)

# plotting --------------------------------------------------------------------
plotd=res1_all %>% 
  # mutate(U95CI=Estimate + 1.96*SE,L95CI=Estimate - 1.96*SE) %>%
  mutate(pheno=str_split(res1_all$roi,"_",simplify=T)[,2]) %>%
  mutate(roi=str_split(res1_all$roi,"_",simplify=T)[,1]) %>%
  filter(E4_status != 'APOE4-by-BMI interaction') %>% 
  mutate(E4_status=factor(E4_status,levels=c("E4-noncarrier","E4-carrier"))) %>%
  mutate(group=factor(group, levels=(bmi.strata.lelves)))
# yrange_forest=0.75*c(-1,1)*max(abs(ranges_95CI))
yrange_forest=0.75*c(-1,1)*2.2#** REMOVE **#

cleanup=theme(panel.grid.major.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                panel.grid.major.y=element_line(colour="grey90",linetype=3),
                panel.background=element_blank(),
                axis.line=element_line(colour="black"),
                panel.border=element_rect(colour="grey", fill=NA, linewidth=0.5))

## sample sizes ---------------------------------------------------------------
plotd %>% ggplot(aes(x=group, y=N, color=E4_status, shape=model)) + geom_point() + 
  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))

## forest plot ----------------------------------------------------------------
plot.alpha=0.7
pos <- position_dodge(width=0)#

# for (model_plot in c("base","full")){
for (model_plot in c("base")){
  for (pheno_plot in c("volume","area","thickness") ){
    plotd_i=plotd %>% filter(model==model_plot, pheno==pheno_plot) %>% arrange(roi,group)
    roi.levels=plotd_i %>% filter(E4_status=="E4-carrier",group==sort(unique(plotd_i$group))[length(sort(unique(plotd_i$group)))]) %>% arrange(Estimate)
    roi.levels=roi.levels$roi
    plotd_i=plotd_i %>% mutate(roi=factor(roi,levels=rev(roi.levels))) %>% arrange(roi)
    p <- plotd_i %>%
      ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,
                 # alpha=group,
                 # group=E4_status,
                 color=E4_status,
                 fill=E4_status,
                 group=group,
                 # color=group,
                 # fill=group,
                 shape=group)) + 
      geom_point(size=1.15, stroke=1, position=pos,alpha=plot.alpha) +
      scale_shape_manual(values=c(1,2,6,12)) + 
      geom_hline(yintercept=0, size=0.4) +
      # geom_errorbar(position=pos, width=0.25,size=0.5,alpha=plot.alpha) + 
      # facet_grid(E4_status~.)+
      geom_path(position=pos,alpha=0.5) +
      geom_hline(yintercept=0, linewidth=0.5)+
      xlab(NULL) + 
      coord_flip(ylim=yrange_forest) +
      scale_y_continuous(labels=scales::label_number(accuracy=0.01)) + 
      ggtitle(paste("BMI vs. cortical",pheno_plot,"associations")) + 
      # facet_grid(.~E4_status)+
      # scale_fill_brewer(palette="Dark2") +
      # scale_color_brewer(palette="Dark2") +
      # facet_grid(.~group)+
      facet_grid(.~E4_status)+
      scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) +
      scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) +
      geom_ribbon(aes(ymin=L95CI,ymax=U95CI), alpha=0.15, position=pos, color=NA) #+       guides(col=guide_legend(reverse=F,title=element_blank()), 
             # shape= guide_legend(reverse=F,title=element_blank()),
             # fill=guide_legend(show=FALSE) ) +
    p + 
      cleanup+ 
      labs(color=NULL,fill=NULL)+
      labs(shape=NULL)+
      theme(legend.position='top',legend.box="vertical", legend.margin=margin())+ 
      # theme(plot.margin=margin(t=0.15, r=1.5, b=0.15, l=0.15, "cm"),
      #       axis.text.x=element_text(angle=-30, hjust=0, vjust=0.5)) + 
      ylab("Estimate in SD units")
    
    # ggsave(file.path(outdir,outdir2,paste("forest_plot_by_BMI_classes_",age_group,"_",model_plot,"_",pheno_plot,"_by_E4status.png",sep='')),
    pfile=paste("forest_plot_by_BMI_classes_",
                # age_group,"_",
                model_plot,"_",
                pheno_plot,"_by_BMIclass",Sys.Date(),".png",sep='')
    if(is.null(outdir2)){
      pfile=file.path(outdir,pfile)
    }else{
      pfile=file.path(outdir,outdir2,pfile)
    }
    ggsave(pfile,width=12*0.8,height=8.5*0.8,units='in',dpi=300)
    # width=8.5,height=11,units='in',dpi=300)
  }
}
