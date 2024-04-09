getN = function(data_name, filter_expression){
  execute_string = paste("analdat <- ",data_name,"%>% filter(",filter_expression,")",sep='')
  eval(parse(text=execute_string))
  analdat
  ret = nrow(analdat)
  ret
}

NE4_carrier <- nrow(d %>% filter(E4_status=="E4-carrier"))
NE4_noncarrier <- nrow(d %>% filter(E4_status=="E4-noncarrier"))
filter_is = 2:35
for(filter_i in filter_is){
  filter_expression_i = df_filter$expression[filter_i]
  analdat = get_subset("d",filter_expression = filter_expression_i)
  NE4_carrier = c(NE4_carrier,nrow(analdat %>% filter(E4_status=="E4-carrier")))
  NE4_noncarrier = c(NE4_noncarrier,nrow(analdat %>% filter(E4_status=="E4-noncarrier")))
  rm(analdat)
}

df_filter = df_filter %>% mutate(N_E4_noncarrier=NE4_noncarrier,N_E4_carrier=NE4_carrier)
cat(df_filter$group_name,sep='\n')
