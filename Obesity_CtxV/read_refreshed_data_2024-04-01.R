library(R.utils)
library(data.table)
library(tidyverse)
library(ukbtools)

datdir = '/Users/jshin/OneDrive - SickKids/ukbb/download/ukb677610'
setwd("/Users/jshin/Documents/ukb_analyses/Obesity_CtxV")

my_ukb_data <- ukb_df("ukb677610", path = "/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb/download/ukb677610")
install.packages('RSQLite')
d = fread('/Users/jshin/Library/CloudStorage/OneDrive-SickKids/ukbb/download/ukb677610/ukb677610.tab.gz')

library(XML)
d_header = readHTMLTable("/Users/jshin/OneDrive - SickKids/ukbb/download/ukb677610/ukb677610.html")
head(d_header[[2]])
my_ukb_key <- ukb_df_field("ukb677610", path = datdir)
names(my_ukb_key)
head(my_ukb_key$col.name)
my_ukb_key$col.name[str_detect(my_ukb_key$col.name,"insula")]
ind = str_detect(my_ukb_key$col.name,"volume_of_");print(sum(ind))
ind = ind & str_detect(my_ukb_key$col.name,"_left_hemisphere");print(sum(ind))
