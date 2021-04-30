## 2 calculate library size; stratify by groups
library(dplyr)
options(stringsAsFactors = F)

setwd('./data/proc/')
cmv = readRDS('cmv_all.rds') # 709 cmv
covid = readRDS('covid_all.rds') # 1477 covid

## library size
cmv.lb = sapply(cmv,function(x) unique(x$productive_templates)) %>% unlist
cmv.dat = data.frame(sample_name = names(cmv.lb),productive_templates = cmv.lb) %>% 
    arrange(productive_templates) %>%
    mutate(group = cut(productive_templates,breaks = quantile(cmv.lb,seq(0,1,0.1)),include.lowest = T,labels = paste0('g',1:10)) %>% as.character)

str(cmv.dat)

covid.lb = sapply(covid,function(x) unique(x$productive_templates)) %>% unlist
covid.dat = data.frame(sample_name = names(covid.lb),productive_templates = covid.lb) %>% 
    arrange(productive_templates) %>%
    mutate(group = cut(productive_templates,breaks = quantile(covid.lb,seq(0,1,0.1)),include.lowest = T,labels = paste0('g',1:10)) %>% as.character)
str(covid.dat)

comb.dat = rbind(cmv.dat,covid.dat) %>% arrange(group,productive_templates)
smallest = comb.dat %>% group_by(group) %>% summarise(smallest = min(productive_templates))
comb.dat$downsample = ifelse(comb.dat$productive_templates %in% smallest$smallest,F,T)
table(comb.dat$group)
table(comb.dat$downsample)
write.csv(comb.dat,'librarysize.csv',row.names = F)
