setwd('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/Grant/NIH_2022Nov_ZP/neuroCHARGE_scripts')

# BMI-ctxVolume NeuroCHARGE project ------------------------------------------#

# tableOne --------------------------------------------------------------------
# https://ehsanx.github.io/intro2R/data-summary-with-tableone.html
filter_i0 = 12; filter_i1 = 17
analdat = c()
	dplyr::select(IID,BMI,age,sex) %>%
	mutate(group = group_name)	
	analdat = rbind(analdat,analdati)
}
#dput(names(table(analdat$group)))

bmi.levels= c("BMI <20kg/m2", 
"20<= BMI <25kg/m2", 
"25<= BMI <30kg/m2", 
"30<= BMI <35kg/m2", 
"35<= BMI <40kg/m2", 
"BMI >=40kg/m2")

no_strata = CreateTableOne(data = analdat %>% mutate(group=factor(group,levels=bmi.levels)),
                         vars = c("age", "sex", "group"), 
                         factorVars = c("sex","group")
                         )

print(no_strata,
      showAllLevels = TRUE,
      nonnormal = "Age"
     )
tab_csv <- print(no_strata,
                 nonnormal = "age",
                 printToggle = FALSE)
write.csv(tab_csv, file = file.path(outdir,"tableOne_summary_age_sex_BMI_group.csv"))


strata <- CreateTableOne(data = analdat,
                         vars = c("age", "sex"), ## Note that BMI-group is not included because we already have strata = Gender
                         factorVars = c("sex"), ## Again, BMI-group is not included because it is in the strata argument
                         strata = "group" ## BMI-group
                         )
print(strata, nonnormal="age",cramVars = 'group')

tab_csv <- print(strata,
                 nonnormal = "age",
                 printToggle = FALSE)
write.csv(tab_csv, file = file.path(outdir,"tableOne_summary_age_sex_by_BMI_group.csv"))