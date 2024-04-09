#*****************************************************************************#
#
# 3c. (BMI >=25 kg/m2) 
# Compare the association results between base vs. full models.
#
BMI_group = 'BMI_gt_25'
#*****************************************************************************#
op <- options(warn=1)

cat("\n3c. Loading association estimates for BMI vs. roi-ctxTH for both the base and full models (BMI >=25 kg/m2).\n")
logfile=file.path(outdir,paste("3c_",BMI_group,"_plot_roiAssociation_results_base_vs_fullModels.log",sep='')); messages=file(logfile, open="wt")
messages=file(logfile, open="wt")
sink(messages, type="message")
sink(messages, type="output")
resfiles = list(base=file.path(outdir,paste("roi_assoc_res_two_methods_baseModel_",BMI_group,".Rdata",sep='')),
                full=file.path(outdir,paste("roi_assoc_res_two_methods_fullModel_",BMI_group,".Rdata",sep='')))
# forest and ggseg plots ------------------------------------------------------
for(mod in c('base','full')){
  load(resfiles[[mod]])
  ## compare two adjustment methods
  beta.p = roi_assoc_res_two_methods %>% 
    ggplot(aes(x=Estimate_reg_out_covs,y=Estimate_adj_covs_in_mod,color=term)) +
    geom_point() + geom_abline(intercept = 0, slope = 1) + 
    theme_bw() +
    theme(legend.position = 'none')
  SE.p = roi_assoc_res_two_methods %>% 
    ggplot(aes(x=SE_reg_out_covs,y=SE_adj_covs_in_mod,color=term)) +
    geom_point() + geom_abline(intercept = 0, slope = 1)+ 
    theme_bw() +
    theme(legend.position = 'none')
  Pval.p = roi_assoc_res_two_methods %>% 
    ggplot(aes(x=-log10(P_reg_out_covs),y=-log10(P_adj_covs_in_mod),color=term)) +
    geom_point() +   theme_bw() +
    geom_abline(intercept = 0, slope = 1)
  p_method_comp = beta.p + SE.p + Pval.p
  
  png(file.path(outdir,paste('roi_assoc_res_two_methods_',model,"_",BMI_group,'.png',sep='')),
      width=8.5,height = 3, units="in", res=300)
  print(p_method_comp +  plot_annotation(
    title = 'Comparing two approaches to adjusting for covariates',
    subtitle = 'Method 1 (x-axis): correct first -> examine associations vs. Method2 (y-axis): correct covariates while fitting.'
  ))
  dev.off()
  
  ##
  roi.dk = subset(data.frame(dk$data), hemi=='left', select=c('region','label'))
  roi.dk$label2 = str_remove(roi.dk$label,"lh_")
  table(roi.dk$label2)
  roi.dk = unique(roi.dk)
  
  for(pheno in c("volume",'area','thickness')){
    roi_ctx_measurements.ordered = res1 %>% filter(E4_status=="E4-carrier", str_detect(roi,pheno)) %>% 
      arrange(Estimate)
    roi_ctx_measurements.ordered = roi_ctx_measurements.ordered$roi
    res1_pheno = res1 %>% 
      filter(str_detect(roi,pheno),E4_status %in% c("E4-carrier","E4-noncarrier")) %>%
      mutate(roi=factor(roi,levels=rev(roi_ctx_measurements.ordered))) %>% 
      mutate(U95CI = Estimate + 1.96*SE,L95CI = Estimate - 1.96*SE) %>%
      arrange(roi)
    
    pos <- position_dodge(width=0.75)#
    p <- res1_pheno %>% 
      ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,color=E4_status,group=E4_status)) + 
      geom_point(size=1, stroke=0.5, position=pos,alpha=0.5) +
      geom_errorbar(position=pos, width=0.25,size=0.5,alpha=0.5) + 
      geom_hline(yintercept=0, size=0.4) +
      geom_path(aes(group=E4_status),position=pos,alpha=0.5) +
      geom_hline(yintercept = 0, linewidth=0.5)+
      xlab(NULL) + 
      coord_flip() +
      ylim((range(res1_pheno$L95CI,res1_pheno$U95CI)))+#only when coord_flip
      ggtitle(paste("BMI vs. cortical",pheno,"associations")) + 
      scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) + 
      guides(col = guide_legend(reverse=T))
    
    pdf(file.path(outdir,paste("forest_plot_assoc_BMI_roi_ctx_",pheno,"_",model,"_",BMI_group,".pdf",sep='')),
        width=6,height=6.5)
    print(p + theme_bw())
    dev.off()
    
    ## ggseg
    plotData_all = c()
    for(e4_status in c("E4-carrier","E4-noncarrier")){
      plotData = merge(roi.dk,
                       res1_pheno %>% 
                         mutate(roi = str_remove(roi,paste("_",pheno,sep=''))) %>%
                         filter(E4_status == e4_status),
                       by.x='label2',by.y='roi',all.x=T)
      plotData$term <- "BMI"
      plotData$method = "m1"
      plotData$E4_status = e4_status
      plotData_all = rbind(plotData_all,plotData)
    }#for(e4_status...)
    plotData_all = plotData_all %>% 
      mutate(E4_status = factor(E4_status,levels=c("E4-noncarrier","E4-carrier")))
    scale.limits = range(plotData_all$Estimate,na.rm=T)
    scale.limits = range(c(scale.limits,-scale.limits))
    my.breaks = seq(round(scale.limits[1],3),round(scale.limits[2],3),length.out=5)
    
    p_ggseg = plotData_all %>%
      group_by(E4_status) %>% 
      ggseg(hemisphere="left",
            mapping=aes(fill = Estimate), color="black", size=0.5,
            atlas = dkextra,
            show.legend = TRUE) +
      scale_fill_brain("dk") +
      facet_grid(cols=vars(E4_status)) +
      labs(title=paste("BMI vs. cortical ",pheno, " under ", model, "\n(",BMI_group,")", sep=""),
           fill="Estimate in SD units") +
      scale_fill_gradient2(low = "darkblue", high="darkred", mid="white", limits=scale.limits, 
                           breaks=my.breaks) + #limits = range of effect sizes, breaks are for the legend
      theme(axis.text = element_blank(), 
            axis.title = element_blank(),
            # strip.text = element_blank(),
            strip.background = element_rect(fill=NA, color=NA),
            legend.position = "bottom",
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(1.5, "cm"),
            legend.text.align = 0.5,
            legend.text = element_text(size=8),
            legend.title = element_blank())#element_text(size=10,vjust=0.75))
    png(file.path(outdir,paste("brain_image_ctx_",pheno,"_associations_wi_BMI_byApoE_",model,"_",BMI_group,".png",sep='')),
        height=4, width=6, units="in", res=300)
    print(p_ggseg + theme(panel.spacing = unit(1, "lines")))
    dev.off()
  }
  rm(p,res1,res2,model,BMI_group)
  }

# starting -------------------------------------------------------------------- 
load(resfiles$base)
res1_base = res1
res2_base = res2 
roi_assoc_res_two_methods_base = roi_assoc_res_two_methods

load(resfiles$full)
res1_full = res1
res2_full = res2 
roi_assoc_res_two_methods_full = roi_assoc_res_two_methods

rm(res1,res2,roi_assoc_res_two_methods)

# forest plot
res1 = c()
res1 = rbind(res1,data.frame(res1_base,model="base"))
res1 = rbind(res1,data.frame(res1_full,model="full"))

for(pheno in c('volume','area','thickness')){
  roi.34.ordered = res1_base %>% filter(str_detect(roi,pheno),E4_status=="E4-carrier") %>% arrange(Estimate)
  roi.34.ordered = roi.34.ordered$roi
  res1_pheno = res1 %>% filter(str_detect(roi,pheno), E4_status != 'APOE4-by-BMI interaction') %>%
    mutate(roi = factor(roi,levels=rev(roi.34.ordered)),
           model = factor(model, levels=c('base','full')),
           E4_status = factor(E4_status, levels=c('E4-noncarrier','E4-carrier')),
           U95CI = Estimate + 1.96*SE,
           L95CI = Estimate - 1.96*SE) %>% 
    arrange(roi)
  alpha.scale = 0.5
  pos <- position_dodge(width=0.75)#
  p <- res1_pheno %>% 
    ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,color=E4_status,shape=model)) + 
    geom_point(size=2, stroke=0.5, position=pos,alpha=alpha.scale) +
    geom_errorbar(position=pos, width=0.25,linewidth=0.5,alpha=alpha.scale) + 
    geom_hline(yintercept=0, size=0.4) +
    geom_path(aes(group=model),position=pos,alpha=alpha.scale) +
    scale_shape_manual(values=c('base'=19,'full'=17)) + 
    scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) +
    geom_hline(yintercept = 0)+
    xlab(NULL) + 
    facet_grid(cols = vars(E4_status)) + 
    coord_flip() +
    ylim((range(res1_pheno$L95CI,res1_pheno$U95CI)))+#only when coord_flip
    ggtitle(paste("BMI vs. cortical",pheno,"associations")) + 
    theme_bw() +
    guides(col = guide_legend(reverse=F),
           linetype = guide_legend(reverse=F),
           shape = guide_legend(reverse=F))
  pdf(file.path(outdir,paste("forest_plot_assoc_BMI_roi_ctx_",pheno,"_base_vs_full_Model_",BMI_group,".pdf",sep='')),
      width=15,height=8)
  print(p)
  dev.off()
}

sink()
closeAllConnections()
print(readLines(logfile))

cat("Finished forest plot of association estimates for BMI vs. roi-ctxTH for both the base and full models (BMI >=25 kg/m2).\n")
options(op)
