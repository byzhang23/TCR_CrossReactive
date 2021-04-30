## extract original count and frequency
library(dplyr)
library(data.table)
library(readxl)
library(tidyr)
options(stringsAsFactors = F)

setwd('./covid_adaptive/')

all = read.csv('./data/Sample_Column_addpseudocount.csv')
input.sample = unique(all$sample)

extract_info <- function(input = 'CCP11',
                      input.dir = './data/manafest_readout/',
                      output.dir = './result/manafest_readout/'){
    
    sub = all %>% filter(sample==input)
    sub$column = gsub('abundance','percent',sub$column)
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
    colnames(d) = gsub(paste0('^',input,'_'),'',colnames(d))
    print(range(d))
    
    # (4) combine
    comb = data.frame(amino_acid = rownames(d),sample = unique(sub$sample),d,check.names = F)
    write.csv(comb,paste0(output.dir,input,'_rawdata.csv'),row.names = F)
    
    comb
}

raw.ls = lapply(input.sample,function(s){
    extract_info(input = s,
                 input.dir = './data/manafest_readout/',
                 output.dir = './result/manafest_readout/')
})

output.dir = './result/manafest_readout/'
ttl.columns = setdiff(gsub(paste(paste0(input.sample,'_'),collapse = '|'),'',gsub('abundance','percent',all$column)) %>% unique %>% sort,'HIV_percent')
raw = do.call(rbind,lapply(raw.ls,function(m){
    
    if(sum(colnames(m)=='HIV_percent')) colnames(m)[colnames(m)=='HIV_percent'] = 'HIV_A_percent'
    miss = setdiff(ttl.columns,colnames(m))
    if(length(miss)>0){
        miss.mat = matrix(NA,ncol = length(miss),nrow = nrow(m))
        colnames(miss.mat) = miss
        m = cbind(m,miss.mat)
    }
    
    m[,c('amino_acid','sample',paste0(rep(c('HIV','COV2S','NL63S','OC43S','HKU1S','229ES'),each = 3),'_',rep(c('A','B','C'),3),'_percent'))]
}))
write.csv(raw,paste0(output.dir,'combine_rawdata.csv'),row.names = F)


