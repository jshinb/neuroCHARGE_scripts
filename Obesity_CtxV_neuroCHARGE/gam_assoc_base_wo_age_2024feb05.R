analdat = analdat %>% mutate(logBMI = log(BMI))
pheno = 'medialorbitofrontal thickness'
library(gratia)
load_mgcv()
{
  model_old_noncarrier = gam(
    medialorbitofrontal_thickness.adj~
      sex+
      s(age,bs='cr',k=4)+ 
      s(BMI,bs='cr',k=4),
    data = analdat %>% 
      mutate(E4_status=factor(E4_status,levels=c("E4-noncarrier","E4-carrier")),
             sex = factor(sex,levels=c("M","F")),
             education = factor(education),
             icv.scaled = scale(ICV)[,1]) %>%
      filter(age>=65,E4_status=="E4-noncarrier"),
  )
  plot(model_old_noncarrier,pages=1)
  fd_old_noncarrier <- derivatives(model_old_noncarrier, term="s(BMI)",partial_match = TRUE)  
  #  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
  #  scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))
  
  p_old_noncarrier = draw(fd_old_noncarrier) + geom_hline(yintercept = 0, linewidth=0.2, col="red", linetype=2) + 
    ggtitle("old E4-noncarrier")
}

{
  model_old_carrier = gam(
    medialorbitofrontal_thickness.adj~
      sex+
      s(age,bs='cr',k=4)+ 
      s(BMI,bs='cr',k=4),
    data = analdat %>% 
      mutate(E4_status=factor(E4_status,levels=c("E4-noncarrier","E4-carrier")),
             sex = factor(sex,levels=c("M","F")),
             education = factor(education),
             icv.scaled = scale(ICV)[,1]) %>%
      filter(age>=65,E4_status=="E4-carrier"),
  )
  plot(model_old_carrier,pages = 1)
  fd_old_carrier <- derivatives(model_old_carrier, term="s(BMI)",partial_match = TRUE)  
  #  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
  #  scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))
  
  p_old_carrier = draw(fd_old_carrier) + geom_hline(yintercept = 0, linewidth=0.2, col="red", linetype=2) + 
    ggtitle("old E4-carrier")
}

{
  model_young_noncarrier = gam(
    medialorbitofrontal_thickness.adj~
      sex+
      s(age,bs='cr',k=4)+ 
      s(BMI,bs='cr',k=4),
    data = analdat %>% 
      mutate(E4_status=factor(E4_status,levels=c("E4-noncarrier","E4-carrier")),
             sex = factor(sex,levels=c("M","F")),
             education = factor(education),
             icv.scaled = scale(ICV)[,1]) %>%
      filter(age<65,E4_status=="E4-noncarrier"),
  )
  plot(model_young_noncarrier,pages = 1)
  fd_young_noncarrier <- derivatives(model_young_noncarrier, term="s(BMI)",partial_match = TRUE)  
  #  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
  #  scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))
  
  p_young_noncarrier = draw(fd_young_noncarrier) + geom_hline(yintercept = 0, linewidth=0.2, col="red", linetype=2) + 
    ggtitle("young E4-noncarrier")
}

{
  model_young_carrier = gam(
    medialorbitofrontal_thickness.adj~
      sex+
      s(age,bs='cr',k=4)+ 
      s(BMI,bs='cr',k=4),
    data = analdat %>% 
      mutate(E4_status=factor(E4_status,levels=c("E4-noncarrier","E4-carrier")),
             sex = factor(sex,levels=c("M","F")),
             education = factor(education),
             icv.scaled = scale(ICV)[,1]) %>%
      filter(age<65,E4_status=="E4-carrier"),
  )
  plot(model_young_carrier,pages = 1)
  fd_young_carrier <- derivatives(model_young_carrier, term="s(BMI)",partial_match = TRUE)  
  #  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
  #  scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))
  
  p_young_carrier = draw(fd_young_carrier) + geom_hline(yintercept = 0, linewidth=0.2, col="red", linetype=2) + 
    ggtitle("young E4-carrier")
}

p1=draw(model_young_noncarrier, select = "BMI", partial_match = T)+ggtitle("young E4-noncarrier")
p2=draw(model_young_carrier, select = "BMI", partial_match = T)+ggtitle("young E4-carrier")
p3=draw(model_old_noncarrier, select = "BMI", partial_match = T)+ggtitle("old E4-noncarrier")
p4=draw(model_old_carrier, select = "BMI", partial_match = T)+ggtitle("old E4-carrier")

p_xlim = range(analdat$BMI,na.rm=T)
p_ylim1= c(-1.75,0.5)#c(-1.25,0.5)
p_ylim2= c(-0.075,0.06)#c(-2.25,2)
# p1+p3+p2+p4 + plot_layout(ncol = 2, nrow=2) & coord_cartesian(xlim=p_xlim,ylim=p_ylim1)&
#   plot_annotation(title = pheno,
#                   theme = theme(plot.title = element_text(size = 12)))
# 
p1=p1+scale_y_continuous(limits=p_ylim1)
p_young_noncarrier=p_young_noncarrier+scale_y_continuous(limits=p_ylim2)
p3=p3+scale_y_continuous(limits=p_ylim1)
p_old_noncarrier=p_old_noncarrier+ scale_y_continuous(limits=p_ylim2)
p2=p2+scale_y_continuous(limits=p_ylim1)
p_young_carrier=p_young_carrier+scale_y_continuous(limits=p_ylim2)
p4=p4+scale_y_continuous(limits=p_ylim1)
p_old_carrier=p_old_carrier+scale_y_continuous(limits=p_ylim2)
p1+p3+p2+p4+
  p_young_noncarrier+
  p_old_noncarrier+
  p_young_carrier+
  p_old_carrier+
plot_layout(ncol = 4, nrow=2) & coord_cartesian(xlim=p_xlim)&
  plot_annotation(title = pheno,
                  theme = theme(plot.title = element_text(size = 12)))

# ggsave(paste(pheno,"_gam_derivatives_wo_age_",Sys.Date(),".png",sep=),
#        width=20*0.7,height=11*0.7,dpi=300)
ggsave(paste(pheno,"_gam_derivatives_wi_age_",Sys.Date(),".png",sep=),
       width=20*0.7,height=11*0.7,dpi=300)
