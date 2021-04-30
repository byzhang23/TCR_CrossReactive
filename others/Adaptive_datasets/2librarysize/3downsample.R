## 3-downsampling
library(dplyr)
options(stringsAsFactors = F)
setwd()

output.dir = './data/proc/downsample/'
if(!dir.exists(output.dir)) dir.create(output.dir)

## read in dataset
meta = read.csv('./data/proc/librarysize.csv') 

downsample <- function(dat,
                       tg.group,
                       seed = 1){
    
    ref = meta %>% filter(group==tg.group)
    
    remain = ref$sample_name[ref$downsample==T]
    lb.size = ref$productive_templates[ref$downsample==F]
    
    sub = dat[names(dat) %in% remain]
    sub.down = lapply(1:length(sub),function(i){
        
        print(i)
        dt = sub[[i]]
        dt$id = 1:nrow(dt)
        set.seed(12345 + seed)
        id = sample(1:nrow(dt),size = lb.size,replace = T,prob = dt$productive_frequency)
        new = table(id) %>% data.frame %>% 
            mutate(id = as.integer(id),productive_templates = lb.size) %>% 
            dplyr::rename(templates = Freq) %>% 
            mutate(productive_frequency = templates/lb.size)
        
        final = inner_join(dt %>% select(id,sample_name,amino_acid),new) %>% select(-id) %>% filter(!is.na(templates))
        
    })
    names(sub.down) = names(sub)
    
    ## final: # smallest (no downsample) + remaining
    if(ref$sample_name[ref$downsample==F] %in% names(dat)){
        res = dat[names(dat)==ref$sample_name[ref$downsample==F]]
        final = append(sub.down,res)
        names(final) = c(names(sub.down),ref$sample_name[ref$downsample==F])
    }else{
            final = sub.down
        }
    return(final)
    
}

## cmv
cmv = readRDS('./data/proc/cmv_all.rds')
covid = readRDS('./data/proc/covid_all.rds')

for(b in 1:20){
    
    print(paste0('###### ',b,' #####'))
    final.cmv = list()
    for(tg.group in unique(meta$group)){
        print(paste0('====== cmv_',tg.group,' ======='))
        final.cmv = append(final.cmv,
                           downsample(cmv,tg.group,b))
    }
    print(length(final.cmv)==709)
    saveRDS(final.cmv,paste0(output.dir,'cmv_',b,'.rds'))
    
    ## covid
    final.covid = list()
    for(tg.group in unique(meta$group)){
        
        print(paste0('====== covid_',tg.group,' ======='))
        final.covid = append(final.covid,
                             downsample(covid,tg.group,b))
    }
    print(length(final.covid)==1477)
    saveRDS(final.covid,paste0(output.dir,'covid_',b,'.rds'))
    
}



