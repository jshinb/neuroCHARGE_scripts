library("gratia")
library("mgcv")
library("ggplot2")
library("dplyr")
library("stringr")
library(patchwork)
# age- and sex-corrected
# analdat = analdat %>% mutate(logBMI = log(BMI))

load_mgcv()
analdat = d %>% mutate(logBMI = log(BMI)) %>% 
  mutate(total_thickness = global_mean_thickness)
rois_org = rois
rois = c('total',rois_org)
for (phenoi in c('area','thickness','volume')){
  # analdat = d %>% mutate(logBMI = log(BMI))
  for(i in 1 ){#length(rois)){
    ylim1_range_roi <- ylim2_range_roi <- c()
    pheno = paste(rois[i],phenoi)
    analdat[['y']] = analdat[[str_replace(pheno,"[ ]","_")]]
    
    #model_old_noncarrier----------------------------------------------------------
    {
      model_old_noncarrier = gam(
        y~
          sex+
          education+s(icv.scaled, bs='cr', k=4) +  
          s(age,bs='cr',k=4)+ 
          s(BMI,bs='cr',k=4),
        data = analdat %>% 
          mutate(E4_status=factor(E4_status,levels=c("E4-noncarrier","E4-carrier")),
                 sex = factor(sex,levels=c("M","F")),
                 education = factor(education),
                 icv.scaled = scale(ICV)[,1]) %>%
          filter(age>=65,E4_status=="E4-noncarrier"),
      )
      gratia::draw(model_old_noncarrier,nrow=1, select='s(BMI)', partial_match = TRUE)
      gratia::draw(model_old_noncarrier,nrow=1, select='s(age)', partial_match = TRUE)
      fd_old_noncarrier <- derivatives(model_old_noncarrier, select="s(BMI)",
                                       partial_match = TRUE)  
      gratia::draw(model_old_noncarrier) + plot_annotation(title=pheno)
      #  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
      #  scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))
      p_old_noncarrier = gratia::draw(fd_old_noncarrier) + geom_hline(yintercept = 0, linewidth=0.2, col="red", linetype=2) + 
        ggtitle(paste("old E4-noncarrier (n=",nobs(model_old_noncarrier),")",sep=''))
      
      summary(model_old_noncarrier)
      ylim1_range = confint(model_old_noncarrier,parm="BMI",type="confidence",partial_match = T)%>% arrange(`.lower_ci`)
      ylim2_range = range(c(fd_old_noncarrier$.lower_ci,fd_old_noncarrier$.upper_ci))
      ylim1_range = range(ylim1_range%>% dplyr::select(matches("_ci")))
      ylim1_range_roi = c(ylim1_range_roi,ylim1_range)
      ylim2_range_roi = c(ylim2_range_roi,ylim2_range)
    }
    
    #model_old_carrier-------------------------------------------------------------
    {
      model_old_carrier = gam(
        y~
          sex+
          education+s(icv.scaled, bs='cr', k=4) +  
          s(age,bs='cr',k=4)+ 
          s(BMI,bs='cr',k=4),
        data = analdat %>% 
          mutate(E4_status=factor(E4_status,levels=c("E4-noncarrier","E4-carrier")),
                 sex = factor(sex,levels=c("M","F")),
                 education = factor(education),
                 icv.scaled = scale(ICV)[,1]) %>%
          filter(age>=65,E4_status=="E4-carrier"),
      )
      draw(model_old_carrier,nrow=1, select='s(BMI)', partial_match = TRUE)
      draw(model_old_carrier,nrow=1, select='s(age)', partial_match = TRUE)
      fd_old_carrier <- derivatives(model_old_carrier, select="s(BMI)",
                                    partial_match = TRUE)  
      draw(model_old_carrier) + plot_annotation(title=pheno)
      #  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
      #  scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))
      
      p_old_carrier = draw(fd_old_carrier) + geom_hline(yintercept = 0, linewidth=0.2, col="red", linetype=2) + 
        ggtitle(paste("old E4-carrier (n=",nobs(model_old_carrier),")",sep=''))
      summary(model_old_carrier)
      ylim1_range = confint(model_old_carrier,parm="BMI",type="confidence",partial_match = T)%>% arrange(`.lower_ci`)
      ylim1_range = range(ylim1_range%>% dplyr::select(matches("_ci")))
      ylim2_range = range(c(fd_old_carrier$.lower_ci,fd_old_carrier$.upper_ci))
      ylim1_range_roi = c(ylim1_range_roi,ylim1_range)
      ylim2_range_roi = c(ylim2_range_roi,ylim2_range)
    }
    
    #model_young_noncarrier--------------------------------------------------------
    {
      model_young_noncarrier = gam(
        y~#.adj~
          sex+
          education+s(icv.scaled, bs='cr', k=4) +  
          s(age,bs='cr',k=4)+ 
          s(BMI,bs='cr',k=4),
        data = analdat %>% 
          mutate(E4_status=factor(E4_status,levels=c("E4-noncarrier","E4-carrier")),
                 sex = factor(sex,levels=c("M","F")),
                 education = factor(education),
                 icv.scaled = scale(ICV)[,1]) %>%
          filter(age<65,E4_status=="E4-noncarrier"),
      )
      draw(model_young_noncarrier) + plot_annotation(title=pheno)
      fd_young_noncarrier <- derivatives(model_young_noncarrier, term="s(BMI)",partial_match = TRUE)  
      #  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
      #  scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))
      
      p_young_noncarrier = draw(fd_young_noncarrier) + geom_hline(yintercept = 0, linewidth=0.2, col="red", linetype=2) + 
        ggtitle(paste("young E4-noncarrier (n=",nobs(model_young_noncarrier),")",sep=''))
      ylim1_range = confint(model_young_noncarrier,parm="BMI",type="confidence",partial_match = T)%>% arrange(`.lower_ci`)
      ylim2_range = range(c(fd_young_noncarrier$.lower_ci,fd_young_noncarrier$.upper_ci))
      ylim1_range = range(ylim1_range%>% dplyr::select(matches("_ci")))
      ylim1_range_roi = c(ylim1_range_roi,ylim1_range)
      ylim2_range_roi = c(ylim2_range_roi,ylim2_range)
    }
    
    #model_young_carrier--------------------------------------------------------
    {
      model_young_carrier = gam(
        y~#.adj~
          sex+
          education+s(icv.scaled, bs='cr', k=4) +  
          s(age,bs='cr',k=4)+ 
          s(BMI,bs='cr',k=4),
        data = analdat %>% 
          mutate(E4_status=factor(E4_status,levels=c("E4-noncarrier","E4-carrier")),
                 sex = factor(sex,levels=c("M","F")),
                 education = factor(education),
                 icv.scaled = scale(ICV)[,1]) %>%
          filter(age<65,E4_status=="E4-carrier"),
      )
      draw(model_young_carrier) + plot_annotation(title=pheno)
      fd_young_carrier <- derivatives(model_young_carrier, term="s(BMI)",partial_match = TRUE)  
      #  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))+
      #  scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3"))
      
      p_young_carrier = draw(fd_young_carrier) + geom_hline(yintercept = 0, linewidth=0.2, col="red", linetype=2) + 
        ggtitle(paste("young E4-carrier (n=",nobs(model_young_carrier),")",sep=''))
      
      ylim1_range = confint(model_young_carrier,parm="BMI",type="confidence",partial_match = T)%>% arrange(`.lower_ci`)
      ylim2_range = range(c(fd_young_carrier$.lower_ci,fd_young_carrier$.upper_ci))
      ylim1_range = range(ylim1_range%>% dplyr::select(matches("_ci")))
      ylim1_range_roi = c(ylim1_range_roi,ylim1_range)
      ylim2_range_roi = c(ylim2_range_roi,ylim2_range)
    }
    
    #plotting --------------------------------------------------------------------
    ## limits
    p_xlim = range(analdat$BMI,na.rm=T)
    p_ylim1 = range(ylim1_range_roi)
    p_ylim2 = range(ylim2_range_roi)
    
    p1=draw(model_young_noncarrier, select = "BMI", partial_match = T)+ggtitle("young E4-noncarrier")
    p2=draw(model_young_carrier, select = "BMI", partial_match = T)+ggtitle("young E4-carrier")
    p3=draw(model_old_noncarrier, select = "BMI", partial_match = T)+ggtitle("old E4-noncarrier")
    p4=draw(model_old_carrier, select = "BMI", partial_match = T)+ggtitle("old E4-carrier")
    
    p1=p1+coord_cartesian(ylim=p_ylim1,xlim=p_xlim)
    p_young_noncarrier=p_young_noncarrier+coord_cartesian(ylim=p_ylim2)
    p3=p3+coord_cartesian(ylim=p_ylim1,xlim=p_xlim)
    p_old_noncarrier=p_old_noncarrier+ coord_cartesian(ylim=p_ylim2,xlim=p_xlim)
    p2=p2+coord_cartesian(ylim=p_ylim1,xlim=p_xlim)
    p_young_carrier=p_young_carrier+coord_cartesian(ylim=p_ylim2,xlim=p_xlim)
    p4=p4+coord_cartesian(ylim=p_ylim1,xlim=p_xlim)
    p_old_carrier=p_old_carrier+coord_cartesian(ylim=p_ylim2,xlim=p_xlim)
    p1+p3+p2+p4+
      p_young_noncarrier+
      p_old_noncarrier+
      p_young_carrier+
      p_old_carrier+
      plot_layout(ncol = 2, nrow=4, byrow=T) & #coord_cartesian(xlim=p_xlim)&
      plot_annotation(title = pheno,
                      theme = theme(plot.title = element_text(size = 12)))
    
    # ggsave(paste(pheno,"_gam_derivatives_wo_age_",Sys.Date(),".png",sep=),
    #        width=20*0.7,height=11*0.7,dpi=300)
    ggsave(paste(pheno,"_gam_derivatives_wi_age_",Sys.Date(),".png",sep=''),
           height=20*0.7,width=11*0.7,dpi=300)
    
    #cleaning ---------------------------------------------------------------------
    analdat = analdat %>% dplyr::select(-y)
  }
}
