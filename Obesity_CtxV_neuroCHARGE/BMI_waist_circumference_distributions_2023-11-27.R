wd = '/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts'
setwd(wd)

pacman::p_load(tidyverse, tidyselect, data.table, readxl, psych, outliers,
               patchwork, ggcorrplot, hrbrthemes, gridExtra, ggseg,
               tableone, FactoMineR, factoextra,
               ppcor, lsmeans, multcomp, ModelMetrics,
               # GGally,caret, 
               Hmisc, pastecs, testit,mgcv,
               futile.logger,tryCatchLog)# for error/warning messages
# UKBB ------------------------------------------------------------------------
input_specification_file = 'cohort_specific_inputs_UKBB.txt'
#
op <- options(nwarnings = 10000, warn=0)
setwd("./") #current directory
source(input_specification_file)

outdir=paste(cohort_name,ancestry,"outputs",sep="_")
dir.create(outdir)
file.copy(input_specification_file,outdir,overwrite = T)

# define functions ------------------------------------------------------------
cat("Prep: sourcing functions\n")
source("0_functions_get_adj_ctxVol_BMI.R")

# 1. data wrangling -----------------------------------------------------------
logfile=file.path(outdir,"1_DataWrangle.log"); 
op = options(warn = 1)
options(keep.source = TRUE)        # source code file name and line number tracking

flog.appender(appender.file(logfile))  # to log into a file instead of console
flog.threshold(INFO)    # TRACE, DEBUG, INFO, WARN, ERROR, FATAL

tryCatchLog(source("1_DataWrangle.R"))
options(op)

# 2. descriptive stats --------------------------------------------------------
source("2_Descriptive_statistics.R")

# 3A. distributions (age-combined) --------------------------------------------
BMI_mean <- d %>% filter(!is.na(E4_status),!is.na(BMI)) %>%
        group_by(E4_status) %>%
        summarise(Mean = mean(BMI))
        
BMI_median <- d %>% filter(!is.na(E4_status),!is.na(BMI)) %>%
        group_by(E4_status) %>%
        summarise(Median = median(BMI))
        
p_BMI = d %>% filter(!is.na(E4_status)) %>% 
ggplot(aes(x=BMI, fill=E4_status, color=E4_status)) +
 geom_histogram(aes(y = ..density..), alpha=0.2, position='identity') +
 scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) +
 scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) + 
 geom_density(alpha = .2) + 
 geom_vline(data = BMI_mean, aes(xintercept = Mean), linetype=c(1,2), 
 color=c("E4-noncarrier"="#7570B3","E4-carrier"="#D95F02"))

d_ukbb = fread("/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb/data/ukbb41448_waist_circumference_weight_height_2024-11-27.tsv") %>%
filter(eid %in% d$IID) #%>%
# dplyr::select(eid,
# waist_circumference_f48_2_0,# waist circumference
# weight_preimaging_f12143_2_0, #weight pre-imaging
# height_f12144_2_0  )# height pre-imaging

p1 = d_ukbb %>% ggplot(aes(x=weight_preimaging_f12143_2_0,y=weight_f21002_2_0)) + geom_point()
p2 = d_ukbb %>% ggplot(aes(x=weight_preimaging_f12143_2_0,y=height_f12144_2_0)) + geom_point()
p3 = d_ukbb %>% ggplot(aes(x=weight_f21002_2_0,y=height_f12144_2_0)) + geom_point()

p1+p2+p3

lm.fit = lm(log(waist_circumference_f48_2_0)~ height_f12144_2_0 + sex, data = left_join(d_ukbb,d %>% dplyr::rename(eid = IID)))

dM = left_join(d_ukbb,d %>% dplyr::rename(eid = IID)) %>% filter(sex=="M")
gam.fit.M = gam(log(waist_circumference_f48_2_0)~ s(height_f12144_2_0), data = dM, na.action=na.exclude)

dF = left_join(d_ukbb,d %>% dplyr::rename(eid = IID)) %>% filter(sex=="F")
gam.fit.F = gam(log(waist_circumference_f48_2_0)~ s(height_f12144_2_0), data = dF, na.action=na.exclude)

dM = dM %>% mutate(waist_adjHeight = resid(gam.fit.M)+coef(gam.fit.M)[1])
dF = dF %>% mutate(waist_adjHeight = resid(gam.fit.F)+coef(gam.fit.F)[1])

d2 = rbind(dM,dF)
log_waist_mean <- d2 %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
        group_by(E4_status) %>%
        summarise(Mean = mean(waist_adjHeight))
        
log_waist_median <- d2 %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
        group_by(E4_status) %>%
        summarise(Median = median(waist_adjHeight))
        
waist_mean <- d2 %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
        group_by(E4_status) %>%
        summarise(Mean = mean(exp(waist_adjHeight)))
        
waist_median <- d2 %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
        group_by(E4_status) %>%
        summarise(Median = median(exp(waist_adjHeight)))

p_waist = d2 %>% filter(!is.na(E4_status)) %>% 
ggplot(aes(x=exp(waist_adjHeight), fill=E4_status, color=E4_status)) +
 geom_histogram(aes(y = ..density..), alpha=0.2, position='identity') +
 scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) +
 scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) + 
 geom_density(alpha = .2) + 
 geom_vline(data = waist_mean, aes(xintercept = Mean), linetype=c(1,2), 
 color=c("E4-noncarrier"="#7570B3","E4-carrier"="#D95F02")) + 
 xlab("waist circumference_adjHeight (cm)") #+ scale_x_continuous(trans='log')

# 3B. distributions (age >=65 years) ------------------------------------------
BMI_mean_older <- d %>% 
  filter(age >=65) %>%
  filter(!is.na(E4_status),!is.na(BMI)) %>%
  group_by(E4_status) %>%
  summarise(Mean = mean(BMI))
# A tibble: 2 × 2
# E4_status      Mean
# <fct>         <dbl>
#   1 E4-noncarrier  26.4
#   2 E4-carrier     26.1

BMI_median_older <- d %>% 
  filter(age >=65) %>%
  filter(!is.na(E4_status),!is.na(BMI)) %>%
  group_by(E4_status) %>%
  summarise(Median = median(BMI))

BMI_median_older
#  A tibble: 2 × 2
#  E4_status     Median
#  <fct>          <dbl>
# 1 E4-noncarrier   25.8
# 2 E4-carrier      25.7

p_BMI_older = d %>% 
  filter(age >=65) %>%
  filter(!is.na(E4_status)) %>% 
  ggplot(aes(x=BMI, fill=E4_status, color=E4_status)) +
  geom_histogram(aes(y = ..density..), alpha=0.2, position='identity') +
  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) +
  scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) + 
  geom_density(alpha = .2) + 
  geom_vline(data = BMI_mean, aes(xintercept = Mean), linetype=c(1,2), 
             color=c("E4-noncarrier"="#7570B3","E4-carrier"="#D95F02"))

dM = left_join(d_ukbb,d %>% dplyr::rename(eid = IID)) %>% filter(sex=="M") #%>% filter(age >=65)
gam.fit.M = gam(log(waist_circumference_f48_2_0)~ s(height_f12144_2_0), data = dM, na.action=na.exclude)

dF = left_join(d_ukbb,d %>% dplyr::rename(eid = IID)) %>% filter(sex=="F")# %>% filter(age >=65)
gam.fit.F = gam(log(waist_circumference_f48_2_0)~ s(height_f12144_2_0), data = dF, na.action=na.exclude)

dM = dM %>% mutate(waist_adjHeight = resid(gam.fit.M)+coef(gam.fit.M)[1])
dF = dF %>% mutate(waist_adjHeight = resid(gam.fit.F)+coef(gam.fit.F)[1])

d2_older = rbind(dM,dF) %>% filter(age>=65)
log_waist_mean_older <- d2_older %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
  group_by(E4_status) %>%
  summarise(Mean = mean(waist_adjHeight))

log_waist_median_older <- d2_older %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
  group_by(E4_status) %>%
  summarise(Median = median(waist_adjHeight))

waist_mean_older <- d2_older %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
  group_by(E4_status) %>%
  summarise(Mean = mean(exp(waist_adjHeight)))
waist_mean_older
# A tibble: 2 × 2
# E4_status      Mean
# <fct>         <dbl>
# 1 E4-noncarrier  89.2
# 2 E4-carrier     88.4

waist_median_older <- d2_older %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
  group_by(E4_status) %>%
  summarise(Median = median(exp(waist_adjHeight)))
waist_median_older
# E4_status     Median
# <fct>          <dbl>
# 1 E4-noncarrier   89.1
# 2 E4-carrier      88.3

p_waist_older = d2_older %>% filter(!is.na(E4_status)) %>% 
  ggplot(aes(x=exp(waist_adjHeight), fill=E4_status, color=E4_status)) +
  geom_histogram(aes(y = ..density..), alpha=0.2, position='identity') +
  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) +
  scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) + 
  geom_density(alpha = .2) + 
  geom_vline(data = waist_mean, aes(xintercept = Mean), linetype=c(1,2), 
             color=c("E4-noncarrier"="#7570B3","E4-carrier"="#D95F02")) + 
  xlab("waist circumference_adjHeight (cm)") #+ scale_x_continuous(trans='log')

# 3C. distributions (age <65 years) --------------------------------------------
BMI_mean_younger <- d %>% 
  filter(age <65) %>%
  filter(!is.na(E4_status),!is.na(BMI)) %>%
  group_by(E4_status) %>%
  summarise(Mean = mean(BMI))
# A tibble: 2 × 2
# E4_status      Mean
# <fct>         <dbl>
#   1 E4-noncarrier  26.7
#   2 E4-carrier     26.6

BMI_median_younger <- d %>% 
  filter(age <65) %>%
  filter(!is.na(E4_status),!is.na(BMI)) %>%
  group_by(E4_status) %>%
  summarise(Median = median(BMI))

BMI_median_younger
#  A tibble: 2 × 2
#  E4_status     Median
#  <fct>          <dbl>
# 1 E4-noncarrier   25.9
# 2 E4-carrier      25.9

p_BMI_younger = d %>% 
  filter(age <65) %>%
  filter(!is.na(E4_status)) %>% 
  ggplot(aes(x=BMI, fill=E4_status, color=E4_status)) +
  geom_histogram(aes(y = ..density..), alpha=0.2, position='identity') +
  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) +
  scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) + 
  geom_density(alpha = .2) + 
  geom_vline(data = BMI_mean, aes(xintercept = Mean), linetype=c(1,2), 
             color=c("E4-noncarrier"="#7570B3","E4-carrier"="#D95F02"))

dM = left_join(d_ukbb,d %>% dplyr::rename(eid = IID)) %>% filter(sex=="M") #%>% filter(age <65)
gam.fit.M = gam(log(waist_circumference_f48_2_0)~ s(height_f12144_2_0), data = dM, na.action=na.exclude)

dF = left_join(d_ukbb,d %>% dplyr::rename(eid = IID)) %>% filter(sex=="F") #%>% filter(age <65)
gam.fit.F = gam(log(waist_circumference_f48_2_0)~ s(height_f12144_2_0), data = dF, na.action=na.exclude)

dM = dM %>% mutate(waist_adjHeight = resid(gam.fit.M)+coef(gam.fit.M)[1])
dF = dF %>% mutate(waist_adjHeight = resid(gam.fit.F)+coef(gam.fit.F)[1])

d2_younger = rbind(dM,dF) %>% filter(age <65)
log_waist_mean_younger <- d2_younger %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
  group_by(E4_status) %>%
  summarise(Mean = mean(waist_adjHeight))

log_waist_median_younger <- d2_younger %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
  group_by(E4_status) %>%
  summarise(Median = median(waist_adjHeight))

waist_mean_younger <- d2_younger %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
  group_by(E4_status) %>%
  summarise(Mean = mean(exp(waist_adjHeight)))
waist_mean_younger
# E4_status      Mean
# <fct>         <dbl>
# 1 E4-noncarrier  87.5
# 2 E4-carrier     87.2
waist_median_younger <- d2_younger %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
  group_by(E4_status) %>%
  summarise(Median = median(exp(waist_adjHeight)))
waist_median_younger
# E4_status     Median
# <fct>          <dbl>
# 1 E4-noncarrier   86.9
# 2 E4-carrier      86.5

p_waist_younger = d2_younger %>% filter(!is.na(E4_status)) %>% 
  ggplot(aes(x=exp(waist_adjHeight), fill=E4_status, color=E4_status)) +
  geom_histogram(aes(y = ..density..), alpha=0.2, position='identity') +
  scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) +
  scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) + 
  geom_density(alpha = .2) + 
  geom_vline(data = waist_mean, aes(xintercept = Mean), linetype=c(1,2), 
             color=c("E4-noncarrier"="#7570B3","E4-carrier"="#D95F02")) + 
  xlab("waist circumference_adjHeight (cm)") #+ scale_x_continuous(trans='log')

xlim_BMI_age_stratified = range(c(d2_older$BMI,d2_younger$BMI),na.rm=T)
xlim_waist_age_stratified = range(c(exp(d2_younger[["waist_adjHeight"]]),exp(d2_older[["waist_adjHeight"]])),na.rm=T)
(p_BMI_younger + coord_cartesian(xlim=xlim_BMI_age_stratified))  + 
  (p_waist_younger + coord_cartesian(xlim=xlim_waist_age_stratified)) + 
  (p_BMI_older+coord_cartesian(xlim=xlim_BMI_age_stratified)) + 
  (p_waist_older + coord_cartesian(xlim=xlim_waist_age_stratified)) + plot_layout(guides="collect", ncol=2)

ggsave("ukbb_bmi_waist_distbn_by_E4status_age_stratified.png",width=11,height=7.5,
       units="in",dpi=300)

# SPS -------------------------------------------------------------------------
input_specification_file = 'cohort_specific_inputs_SPS.txt'
source(input_specification_file)

outdir=paste(cohort_name,ancestry,"outputs",sep="_")
dir.create(outdir)
file.copy(input_specification_file,outdir,overwrite = T)

# define functions ------------------------------------------------------------
replace.neg.wi.na = function(x){
	ind = !is.na(x)
	ind = ind & x<0
	x[ind] = NA
	x
} 
cat("Prep: sourcing functions\n")
source("0_functions_get_adj_ctxVol_BMI.R")

# 1. data wrangling -----------------------------------------------------------
logfile=file.path(outdir,"1_DataWrangle.log"); 
op = options(warn = 1)
options(keep.source = TRUE)        # source code file name and line number tracking

flog.appender(appender.file(logfile))  # to log into a file instead of console
flog.threshold(INFO)    # TRACE, DEBUG, INFO, WARN, ERROR, FATAL

tryCatchLog(source("1_DataWrangle.R"))
options(op)

# 2. descriptive stats --------------------------------------------------------
source("2_Descriptive_statistics.R")
nbins=20
BMI_mean_SPS <- d %>% filter(!is.na(E4_status),!is.na(BMI)) %>%
        group_by(E4_status) %>%
        summarise(Mean = mean(BMI))
        
BMI_median_SPS <- d %>% filter(!is.na(E4_status),!is.na(BMI)) %>%
        group_by(E4_status) %>%
        summarise(Median = median(BMI))
        
p_BMI_SPS = d %>% filter(!is.na(E4_status)) %>% 
ggplot(aes(x=BMI, fill=E4_status, color=E4_status)) +
 geom_histogram(aes(y = ..density..), alpha=0.2, position='identity',bins=nbins) +
 scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) +
 scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) + 
 geom_density(alpha = .2) + 
 geom_vline(data = BMI_mean_SPS, aes(xintercept = Mean), linetype=c(1,2), 
 color=c("E4-noncarrier"="#7570B3","E4-carrier"="#D95F02"))

d_SPS = fread("/Users/jshin/Library/CloudStorage/OneDrive-SickKids/SYS_Database/jean_shin_dataanalysis_completespssubjects_final_20231127_154229.txt")
d_SPS = data.frame(d_SPS)
d_SPS[,-c(1:3)] = apply(d_SPS[,-c(1:3)],2,as.numeric)
d_SPS[,-c(1:3)] = apply(d_SPS[,-c(1:3)],2,replace.neg.wi.na)
psych::describe(d_SPS)
d_SPS = d_SPS%>% filter(uniqueID %in% d$IID) 
d_SPS = d_SPS %>% dplyr::rename(IID=uniqueID) %>% inner_join(d)
names(d_SPS) = str_remove(names(d_SPS),"SPS.N06f_body_anthropometry.")

dM = d_SPS %>% filter(sex=="M");print(dim(dM))#307
gam.fit.M = gam(log(waist)~ s(height), data = dM, na.action=na.exclude)

dF = d_SPS %>% filter(sex=="F");print(dim(dF))#363
gam.fit.F = gam(log(waist)~ s(height), data = dF, na.action=na.exclude)

dM = dM %>% mutate(waist_adjHeight = resid(gam.fit.M)+coef(gam.fit.M)[1])
dF = dF %>% mutate(waist_adjHeight = resid(gam.fit.F)+coef(gam.fit.F)[1])

d2_SPS = rbind(dM,dF)
log_waist_mean <- d2_SPS %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
        group_by(E4_status) %>%
        summarise(Mean = mean(waist_adjHeight))
        
log_waist_median <- d2_SPS %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
        group_by(E4_status) %>%
        summarise(Median = median(waist_adjHeight))
        
waist_mean_SPS <- d2_SPS %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
        group_by(E4_status) %>%
        summarise(Mean = mean(exp(waist_adjHeight)))
        
waist_median_SPS <- d2_SPS %>% filter(!is.na(E4_status),!is.na(waist_adjHeight)) %>%
        group_by(E4_status) %>%
        summarise(Median = median(exp(waist_adjHeight)))

p_waist_SPS = d2_SPS %>% filter(!is.na(E4_status)) %>% 
ggplot(aes(x=exp(waist_adjHeight), fill=E4_status, color=E4_status)) +
 geom_histogram(aes(y = ..density..), alpha=0.2, position='identity',bins=nbins) +
 scale_color_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) +
 scale_fill_manual(values=c("E4-carrier"="#D95F02","E4-noncarrier"="#7570B3")) + 
 geom_density(alpha = .2) + 
 geom_vline(data = waist_mean_SPS, aes(xintercept = Mean), linetype=c(1,2), 
 color=c("E4-noncarrier"="#7570B3","E4-carrier"="#D95F02")) + 
 xlab("waist circumference_adjHeight (cm)") #+ scale_x_continuous(trans='log')

xlim_BMI = range(c(d2$BMI,d2_SPS$BMI),na.rm=T)
xlim_waist = range(c(exp(d2[["waist_adjHeight"]]),exp(d2_SPS[["waist_adjHeight"]])),na.rm=T)
(p_BMI_SPS + coord_cartesian(xlim=xlim_BMI))  + 
(p_waist_SPS + coord_cartesian(xlim=xlim_waist)) + 
(p_BMI+coord_cartesian(xlim=xlim_BMI)) + 
(p_waist + coord_cartesian(xlim=xlim_waist)) + plot_layout(guides="collect", ncol=2)

ggsave("sps_ukbb_bmi_waist_distbn_by_E4status.png",width=11,height=7.5,
units="in",dpi=300)
waist_mean
