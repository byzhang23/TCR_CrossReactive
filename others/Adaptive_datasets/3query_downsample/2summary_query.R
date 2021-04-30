## 3. Aggregate at AA level and summarize info
library(dplyr)
library(data.table)
library(tidyr)
options(stringsAsFactors = F)
setwd('')

for(filename in c('1CrossReactive','2COV2S_Specific','3COV2N_Specific')){ #,'4NL63S_Specific'
    #filename = commandArgs(trailingOnly = T)[1] %>% as.character
    path = paste0('./result/query_downsample/',filename,'/')
    
    full = read.csv(paste0('./data/query/',filename,'.csv')) %>% select(amino_acid) %>% unique
    
    summary = lapply(1:100,function(f){
        
        dt = readRDS(paste0(path,'cmv_',f,'.rds'))
        cmv = dt %>% group_by(amino_acid) %>% summarise(nsample = n_distinct(sample_name)) %>% mutate(dataset = 'cmv')
        
        dt2 = readRDS(paste0(path,'covid_',f,'.rds'))
        covid = dt2 %>% group_by(amino_acid) %>% summarise(nsample = n_distinct(sample_name)) %>% mutate(dataset = 'covid')
        
        all = rbind(cmv,covid) %>% spread(dataset,nsample)
        res = left_join(full,all)
        
    })
    
    names(summary) = paste0('dataet_',1:100)
    saveRDS(summary,paste0('./result/query_downsample/',filename,'_summary.rds'))
}
