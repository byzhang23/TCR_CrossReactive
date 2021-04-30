library(dplyr)
library(data.table)
library(parallel)
options(stringsAsFactors = F)

setwd()

#############################
## ImmuneCODE-Review-002
#############################
## limit to amino_acid !='' & amino_acid !='na'
files = list.files('./data/ImmuneCODE-Review-002/','tsv') #
output.dir = './data/ImmuneCODE-Review-002/recal_freq/' #
ref = fread('./data/SampleOverview_covid.tsv') %>% select(sample_name,productive_templates)

fail = c()
dt.ls = NULL
for(i in 1:length(files)){
    f = files[i]
    print(f)
    dt = fread(paste0('./data/ImmuneCODE-Review-002/',f)) %>%  #
        mutate(sample_name = gsub('\\.tsv','',f)) %>%
        select(sample_name,amino_acid,templates,productive_frequency) %>% 
        filter(amino_acid !='na' & amino_acid !='' & productive_frequency!='na') %>% 
        mutate(productive_frequency = as.numeric(productive_frequency)) %>% 
        group_by(sample_name,amino_acid) %>% 
        summarise(templates = sum(templates),
                  productive_frequency = sum(productive_frequency)) %>% 
        arrange(-templates) %>% 
        mutate(productive_templates = templates[1]/productive_frequency[1])
    saveRDS(dt,paste0(output.dir,gsub('tsv','rds',f)))
    
    # if(!identical(ref$productive_templates[ref$sample_name==unique(dt$sample_name)] %>% as.numeric,
    #               unique(dt$productive_templates) %>% as.numeric)) fail = c(fail,f)
    # 
    # dt.ls[[i]] = dt
    
}
names(dt.ls) = sapply(dt.ls,function(x) unique(x$sample_name)) %>% unlist

# dt.ls = NULL
# for(i in 1:length(files)){
# 
#     f = files[i]
#     print(f)
#     dt.ls[[i]] = readRDS(paste0(output.dir,gsub('tsv','rds',f)))
# }
# names(dt.ls) = sapply(dt.ls,function(x) unique(x$sample_name)) %>% unlist

## check reported library size
calculate.lb = sapply(dt.ls,function(x) unique(x$productive_templates)) %>% unlist
reported.lb = ref$productive_templates;names(reported.lb) = ref$sample_name
match.id = match(names(calculate.lb),names(reported.lb))
reported.lb = reported.lb[match.id]
dat.lb = data.frame(calculate = calculate.lb %>% as.numeric,reported = reported.lb %>% as.numeric) %>% mutate(diff = calculate-reported)

## small lb size
small.lb = which(calculate.lb<10^4)
if(length(small.lb)>0) dt.ls = dt.ls[-small.lb] 
length(dt.ls) %>% print # remaining: 1405 samples
saveRDS(dt.ls,'./data/proc/covid_1.rds')


#############################
## ImmuneCODE-Review-002.2
#############################
## limit to amino_acid !='' & amino_acid !='na'
files = list.files('./data/ImmuneCODE-Review-002.2/','tsv') #
output.dir = './data/ImmuneCODE-Review-002.2/recal_freq/' #
ref = fread('./data/SampleOverview_covid.tsv') %>% select(sample_name,productive_templates)

ex = gsub('rds','tsv',list.files(output.dir,'rds'))
files = setdiff(files,ex)

#fail = c()
#dt.ls = NULL
for(i in 1:length(files)){
    f = files[i]
    print(f)
    dt = fread(paste0('./data/ImmuneCODE-Review-002.2/',f)) %>%  #
        mutate(sample_name = gsub('\\.tsv','',f)) %>%
        select(sample_name,amino_acid,templates,productive_frequency) %>% 
        filter(amino_acid !='na' & amino_acid !='' & productive_frequency!='na') %>% 
        mutate(productive_frequency = as.numeric(productive_frequency)) %>% 
        group_by(sample_name,amino_acid) %>% 
        summarise(templates = sum(templates),
                  productive_frequency = sum(productive_frequency)) %>% 
        arrange(-templates) %>% 
        mutate(productive_templates = templates[1]/productive_frequency[1])
    saveRDS(dt,paste0(output.dir,gsub('tsv','rds',f)))
    
    # if(!identical(ref$productive_templates[ref$sample_name==unique(dt$sample_name)] %>% as.numeric,
    #               unique(dt$productive_templates) %>% as.numeric)) fail = c(fail,f)
    # 
    # dt.ls[[i]] = dt
    
}
names(dt.ls) = sapply(dt.ls,function(x) unique(x$sample_name)) %>% unlist

## check reported library size
calculate.lb = sapply(dt.ls,function(x) unique(x$productive_templates)) %>% unlist
reported.lb = ref$productive_templates;names(reported.lb) = ref$sample_name
match.id = match(names(calculate.lb),names(reported.lb))
reported.lb = reported.lb[match.id]
dat.lb = data.frame(calculate = calculate.lb %>% as.numeric,reported = reported.lb %>% as.numeric) %>% mutate(diff = calculate-reported)


## small lb size
small.lb = which(calculate.lb<10^4)
if(length(small.lb)>0) dt.ls = dt.ls[-small.lb] 

## remove 2 samples with no meta
rm2.id = which(names(dt.ls) %in% c('ADIRP0000257_20200526_Frblood_Repertorie_TCRB','ADIRP0002657_20200618_Frblood_Repertorie_TCRB'))
dt.ls = dt.ls[-rm2.id]
length(dt.ls) %>% print # remaining: 72 samples
saveRDS(dt.ls,'./data/proc/covid_2.rds')
