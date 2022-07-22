#*****************************************************************************#
#
# 3. compare the association results between base vs. full models
#
#*****************************************************************************#
cat("\n3c. Loading association estimates for WMH vs. roi-ctxTH for both the base and full models.\n")
# starting -------------------------------------------------------------------- 
load(file.path(outdir,"roi_assoc_res_two_methods_baseModel.Rdata"))
res1_base = res1
res2_base = res2 
roi_assoc_res_two_methods_base = roi_assoc_res_two_methods

load(file.path(outdir,"roi_assoc_res_two_methods_fullModel.Rdata"))
res1_full = res1
res2_full = res2 
roi_assoc_res_two_methods_full = roi_assoc_res_two_methods

rm(res1,res2,roi_assoc_res_two_methods)

# forest plot
roi.34.ordered = subset(res1_base,sex=="sex-combined") %>% arrange(Estimate)
roi.34.ordered = roi.34.ordered$roi
res1 = c()
res1 = rbind(res1,data.frame(res1_base,model="base"))
res1 = rbind(res1,data.frame(res1_full,model="full"))

res1 = res1 %>% 
  mutate(roi = factor(roi,levels=rev(roi.34.ordered)),
         model = factor(model, levels=rev(c('base','full'))),
         sex = factor(sex, levels=c('sex-combined','M','F')),
         U95CI = Estimate + 1.96*SE,
         L95CI = Estimate - 1.96*SE) %>% 
  arrange(roi)
alpha.scale = 0.9
pos <- position_dodge(width=0.75)#
p <- subset(res1,sex!="sex-by-WMH interaction") %>% 
  ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,linetype=model,color=sex, shape=model)) + 
  geom_point(size=2, stroke=0.5, position=pos,alpha=alpha.scale) +
  geom_errorbar(position=pos, width=0.25,size=0.5,alpha=alpha.scale) + 
  geom_hline(yintercept=0, size=0.4) +
  geom_path(aes(group=model),position=pos,alpha=alpha.scale) +
  scale_linetype_manual(values=c('base'="dashed",'full'='solid')) + 
  scale_shape_manual(values=c('base'=19,'full'=17)) + 
  scale_color_manual(values=c("sex-combined"="#619CFF","M"="#00BA38","F"="#F8766D")) +
  geom_hline(yintercept = 0)+
  xlab(NULL) + 
  facet_grid(cols = vars(sex)) + 
  coord_flip() +
  ylim((range(res1$L95CI,res1$U95CI)))+#only when coord_flip
  ggtitle("WMH vs. cortical thickness associations") + 
  theme_bw() +
  guides(col = guide_legend(reverse=F),
         linetype = guide_legend(reverse=F),
         shape = guide_legend(reverse=F))
pdf(file.path(outdir,"forest_plot_assoc_WMH_roi_ctxTH_base_fullModel.pdf"),
    width=15,height=8)
print(p)
dev.off()
cat("Finished forest plot of association estimates for WMH vs. roi-ctxTH for both the base and full models.\n")
