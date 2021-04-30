## Compare sequence homology

## heatmap for Fig2-boxplots
library(dplyr)
library(data.table)
library(tidyr)
library(pheatmap)
options(stringsAsFactors = F)

setwd('./covid_adaptive/')
output.dir = './revision/'

raw.dir = './data/HIV/'
files = list.files(raw.dir)
dat.ls = lapply(files,function(f){
    
    if(f=="PC5_HIV_A.tsv"){
        
        dt = fread(paste0(raw.dir,f))[,c(4,2)]
        colnames(dt)[1] = 'freq'
        
        dt = dt %>% filter(aminoAcid !='') %>% 
            group_by(aminoAcid) %>% 
            summarise(freq = sum(freq)) %>% 
            arrange(-freq) %>% 
            mutate(sample = sub('\\_.*','',f))
        
    }else{
        dt = fread(paste0(raw.dir,f)) %>% 
            select(freq,cdr3aa) %>% 
            filter(cdr3aa !='') %>% 
            group_by(cdr3aa) %>% 
            summarise(freq = sum(freq)) %>% 
            arrange(-freq) %>% 
            rename(aminoAcid = cdr3aa) %>% 
            mutate(sample = sub('\\_.*','',f))
    }
    
    dt
})
names(dat.ls) = sub('\\_.*','',files)
saveRDS(dat.ls,'./revision/HIV_A.rds')

## cross reactive query
cr = read.csv('./data/query_revision/1CrossReactive.csv')
info = cr %>% group_by(sample) %>% summarise(n = n_distinct(amino_acid))

sampling <- function(dat.ls,
                     info,
                     seed = 1){
    
    
    samp.ls = lapply(dat.ls,function(dt){
        
        n = info$n[info$sample==dt$sample %>% unique]
        set.seed(seed)
        sub = dt[sample(1:nrow(dt),n,replace = F),]
        
        sub %>% select(-freq) %>% rename(amino_acid = aminoAcid) %>% mutate(type = 'HIV') %>% select(type,sample,amino_acid)
    })
    
    do.call(rbind,samp.ls)
}

dat.random = sampling(dat.ls,info,seed = 1)
write.csv(dat.random,'./revision/HIV_seed1.csv',row.names = F)


