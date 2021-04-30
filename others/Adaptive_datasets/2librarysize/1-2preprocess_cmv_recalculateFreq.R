library(dplyr)
library(data.table)
library(parallel)
options(stringsAsFactors = F)

setwd()

# ## 1. cmv dataset - remove counting method ='V1' (76 deleted)
# overview = fread('./data/SampleOverview_cmv.tsv')
# delete = overview$sample_name[overview$counting_method=='v1'] # remove 76 samples; 710 remaining
# for(d in delete){
#     
#     print(d)
#     file.remove(paste0('./data/cmv/',d,'.tsv'))
# }


## 2. For reminaing 710 datasets - recalculate frequency
## limit to amino_acid !='', productive_frequency = templates/productive_templates
files = list.files('./data/cmv/','tsv')
output.dir = './data/cmv/recal_freq/'

dt.ls = mclapply(files,function(f){
    
    print(f)
    dt = fread(paste0('./data/cmv/',f)) %>% 
        select(sample_name,amino_acid,templates,productive_templates,frame_type) %>% 
        filter(amino_acid !='' & frame_type=='In' & templates>0) %>% 
        mutate(productive_templates = sum(templates)) %>% 
        group_by(sample_name,amino_acid,productive_templates) %>% 
        summarise(templates = sum(templates)) %>% 
        arrange(-templates) %>% 
        mutate(productive_frequency = templates/productive_templates)
    saveRDS(dt,paste0(output.dir,gsub('tsv','rds',f)))
    
    return(dt)
},mc.cores = detectCores())
names(dt.ls) = sapply(dt.ls,function(x) unique(x$sample_name)) %>% unlist


# dt.ls = NULL
# for(i in 1:length(files)){
# 
#     f = files[i]
#     print(f)
#     dt.ls[[i]] = readRDS(paste0(output.dir,gsub('tsv','rds',f)))
# }
# names(dt.ls) = sapply(dt.ls,function(x) unique(x$sample_name)) %>% unlist

# check library size
lb.size = sapply(dt.ls,function(x) unique(x$productive_templates)) %>% unlist
range(lb.size)
small.lb = which(lb.size<10^4)
if(length(small.lb)>0) dt.ls = dt.ls[-small.lb] 
length(dt.ls) %>% print # remaining: 709 samples
saveRDS(dt.ls,'./data/proc/cmv_all.rds')

