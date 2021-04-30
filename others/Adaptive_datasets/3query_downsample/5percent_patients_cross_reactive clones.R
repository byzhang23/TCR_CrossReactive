# % patients have at least 1 cross reactive clones
library(dplyr)
setwd('./covid_adaptive/')


dir = '~/Downloads/1CrossReactive/'
cmv.files = paste0('cmv_',1:100,'.rds')
covid.files = paste0('covid_',1:100,'.rds')

summary.cmv = lapply(cmv.files,function(f){
    
    dat = readRDS(paste0(dir,f))
    dat %>% group_by(sample_name) %>% summarise(nAA = n_distinct(amino_acid))
})

summary.covid = lapply(covid.files,function(f){
    print(f)
    dat = readRDS(paste0(dir,f))
    dat %>% group_by(sample_name) %>% summarise(nAA = n_distinct(amino_acid))
})


total.cmv = 709
total.covid = 1477
##### at least 1 cross reactive

counts_cr <- function(data,
                      k){
    all = sapply(data,function(d) nrow(d %>% filter(nAA==k)))
    median(all)
}

cmv.ct = c()
covid.ct = c()
for(k in c(1:8)){
    
    cmv.ct = c(cmv.ct,counts_cr(summary.cmv,k))
    covid.ct = c(covid.ct,counts_cr(summary.covid,k))
    
}
summary = data.frame(k = 0:8,pre_covid = c(total.cmv-sum(cmv.ct),cmv.ct),covid = c(total.covid-sum(covid.ct),covid.ct))
summary$pre_covid_percent = summary$pre_covid/total.cmv
summary$covid_percent = summary$covid/total.covid
write.csv(summary,'./revision/Fig8_piechart.csv',row.names = F)


# cmv1 = sapply(summary.cmv,function(d) nrow(d %>% filter(nAA>=1)))
# median(cmv1) # 455
# median(cmv1)/total.cmv # 64.18%
# 
# covid1 = sapply(summary.covid,function(d) nrow(d %>% filter(nAA>=1)))
# median(covid1) # 805
# median(covid1)/total.covid # 54.5%
# 
# ##### >3 cross reactive
# cmv3 = sapply(summary.cmv,function(d) nrow(d %>% filter(nAA>3)))
# median(cmv3) # 15
# median(cmv3)/total.cmv # 2.12%
# 
# covid3 = sapply(summary.covid,function(d) nrow(d %>% filter(nAA>3)))
# median(covid3) # 35
# median(covid3)/total.covid # 2.37%
