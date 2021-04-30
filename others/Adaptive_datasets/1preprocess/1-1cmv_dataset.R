## 1. combine data into rds
library(data.table)
library(dplyr)
library(parallel)
options(stringsAsFactors = F)

argument = commandArgs(trailingOnly = T)
k = as.integer(as.character(argument[1]))

setwd('')
raw.dir = './data/cmv/'
output.dir = './data/proc/'
if(!dir.exists(output.dir)) dir.create(output.dir)

files = list.files(raw.dir,'.tsv')

## merge every 100 into 1 rds ==========CHECK
start = seq(1,length(files),by=100)
end = ifelse(start+99>length(files),length(files),start+100)

## output
res = do.call(rbind,mclapply(start[k]:end[k],function(s){
    
    print(s)
    if(files[s] !='INCOV067-AC-3_TCRB.tsv'){
        dat = fread(paste0(raw.dir,files[s]),showProgress=F)
        #dat$sample_name = gsub('.tsv','',files[s])
    }else{
        dat = NULL
    }
    dat
    
},mc.cores = 20))
saveRDS(res,paste0(output.dir,'cmv_',k,'.rds'))
