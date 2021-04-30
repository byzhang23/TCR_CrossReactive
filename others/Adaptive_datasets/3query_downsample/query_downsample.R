## query downsample

library(dplyr)
library(data.table)
library(parallel)
options(stringsAsFactors = F)
setwd()

## query file
f = commandArgs(trailingOnly = T)[1] %>% as.character

## merge
dat = readRDS(paste0('./data/proc/downsample/',f))
for(filename in c('1CrossReactive','2COV2S_Specific','3COV2N_Specific')){
    
    print(filename)
    output.dir = paste0('./result/query_downsample/',filename,'/')
    if(!dir.exists(output.dir)) dir.create(output.dir)
    
    query = data.frame(amino_acid = read.csv(paste0('./data/query/',filename,'.csv'))[,3] %>% unique) #type,sample,amino_acid

    summary = do.call(rbind,mclapply(1:length(dat),function(i){
        
        # print(i)
        dt = dat[[i]]
        res = inner_join(query,dt) %>% select(sample_name,
                                               #rearrangement,
                                               amino_acid,
                                               templates,
                                               productive_frequency)
        
        res
        
    },mc.cores = detectCores())) %>% filter(templates>0) %>% mutate(filename = gsub('\\.rds','',f),
                                                                    dataset = sub('\\_.*','',filename),
                                                                    downsample.id = sub('.*\\_','',filename))
    
    saveRDS(summary,paste0(output.dir,f))
    
}

