## add pseudocount
library(dplyr)
library(data.table)
library(readxl)
library(tidyr)
options(stringsAsFactors = F)

setwd('./covid_adaptive/')

addpseudo <- function(input = 'CCP11',
                      input.dir = './data/manafest_readout/',
                      output.dir = './result/manafest_readout/',
                      pseudocount = 1){
    
    sub = all %>% filter(sample==input)
    old = read_excel(paste0(input.dir,input,' FEST Readout.xlsx'),sheet = 'ref_comparison_only')
    colnames(old)[1]='amino_acid'
    
# for specific columns
    dat = old %>% select(any_of(c('amino_acid',sub$column)))
    print(dim(dat))
    print(ncol(dat)-length(sub$column)-1==0)
    
    # (1) remove all 0s
    d = dat %>% select(any_of(sub$column)) %>% as.matrix
    rownames(d) = dat$amino_acid
    d = d[rowSums(d)>0,]
    new = d + pseudocount
    
    # (2) frequency
    ttl = colSums(new)
    freq = t(t(new)/ttl) # cannot directly divide by colsum (something wrong!)
    
    # (3) average
    freq2 = as.data.frame(freq)
    freq2$amino_acid = rownames(freq2)
    long = inner_join(freq2 %>% gather(column,frequency,-amino_acid),sub)
    average = long %>% 
        group_by(amino_acid,type) %>% 
        mutate(type = paste0('mean.',type)) %>% 
        summarise(mean.freq = mean(frequency)) %>% 
        spread(type,mean.freq)
    
    # (4) ratio
    if(sum(colnames(average) %in% 'mean.229ES')>0){average$ratio.229ES = average$mean.229ES/average$mean.HIV}
    if(sum(colnames(average) %in% 'mean.COV2S')>0){average$ratio.COV2S = average$mean.COV2S/average$mean.HIV}
    if(sum(colnames(average) %in% 'mean.HKU1S')>0){average$ratio.HKU1S = average$mean.HKU1S/average$mean.HIV}
    if(sum(colnames(average) %in% 'mean.NL63S')>0){average$ratio.NL63S = average$mean.NL63S/average$mean.HIV}
    if(sum(colnames(average) %in% 'mean.OC43S')>0){average$ratio.OC43S = average$mean.OC43S/average$mean.HIV}
    
    
    # (4) combine
    colnames(freq) = gsub('abundance','frquency',colnames(freq))
    part1 = cbind(new,freq) %>% data.frame %>% mutate(sample = input,amino_acid = rownames(new))
    comb = inner_join(average,part1)
    write.csv(comb,paste0(output.dir,input,'_addpseudo.csv'),row.names = F)
    
    comb
}


## CCPs
all = read.csv('./data/Sample_Column_addpseudocount.csv')
input.sample = unique(all$sample)
final.ls = lapply(input.sample,function(s){
    addpseudo(input = s,
              input.dir = './data/manafest_readout/',
              output.dir = './result/manafest_readout/',
              pseudocount = 1)
})

output.dir = './result/manafest_readout/'
final = do.call(rbind,lapply(final.ls,function(dt){
    
    m = dt[,c(1,grep('^ratio|sample',colnames(dt)))] %>% data.frame
    miss = setdiff(paste0('ratio.',setdiff(unique(all$type),'HIV')),colnames(m))
    if(length(miss)>0){
        miss.mat = matrix(NA,ncol = length(miss),nrow = nrow(m))
        colnames(miss.mat) = miss
        m = cbind(m,miss.mat)
    }
    
    m[,c('amino_acid','sample',"ratio.COV2S", "ratio.NL63S", "ratio.229ES", "ratio.HKU1S", "ratio.OC43S")]
}))
write.csv(final,paste0(output.dir,'combine_addpseudo.csv'),row.names = F)

## PCs
all = read.csv('./data/manafest_readout/PC/Sample_Column_addpseudocount.csv')
input.pc = unique(all$sample)
pc.ls = lapply(input.pc,function(s){
    addpseudo(input = s,
              input.dir = './data/manafest_readout/PC/',
              output.dir = './result/manafest_readout/PC/',
              pseudocount = 1)
})

pc.output.dir = './result/manafest_readout/PC/'
pc = do.call(rbind,lapply(pc.ls,function(dt){
    
    m = dt[,c(1,grep('^ratio|sample',colnames(dt)))] %>% data.frame
    miss = setdiff(paste0('ratio.',setdiff(unique(all$type),'HIV')),colnames(m))
    if(length(miss)>0){
        miss.mat = matrix(NA,ncol = length(miss),nrow = nrow(m))
        colnames(miss.mat) = miss
        m = cbind(m,miss.mat)
    }
    
    m[,c('amino_acid','sample',"ratio.COV2S", "ratio.NL63S", "ratio.229ES", "ratio.HKU1S", "ratio.OC43S")]
}))
write.csv(pc,paste0(pc.output.dir,'combine_addpseudo.csv'),row.names = F)
