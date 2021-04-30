## check library size
library(dplyr)
library(parallel)
options(stringsAsFactors = F)
setwd()

dat.dir = './data/proc/'
output.dir = './result/revision/1librarysize/'
if(!dir.exists(output.dir)) dir.create(output.dir)
# 
# for(initial in c('combine')){
#     
#     print('#######')
#     print(initial)
#     print('#######')
#     
#     res.ls = mclapply(1:8,function(i){
#         
#         print('=====')
#         print(i)
#         print('=====')
#         file = paste0(initial,'_',i,'.rds') # 'combine_','cmv_'
#         r = readRDS(paste0(dat.dir,file))
#         dim(r) %>% print
#         
#         res = r %>% select(templates,productive_frequency,sample_name) %>% 
#             mutate(productive_frequency = as.numeric(productive_frequency),
#                    templates = as.numeric(templates)) %>% 
#             filter(!is.na(productive_frequency)) %>% 
#             group_by(sample_name) %>% 
#             dplyr::summarise(library.size = sum(templates))
#         
#         saveRDS(res,paste0(output.dir,file))
#         
#         return(res)
#     },mc.cores = detectCores())
#     
#     concat = do.call(rbind,res.ls)
#     saveRDS(concat,paste0(output.dir,initial,'.rds'))
# }


## cmv
for(initial in c('cmv')){
    
    print('#######')
    print(initial)
    print('#######')
    
    res.ls = mclapply(1:8,function(i){
        
        print('=====')
        print(i)
        print('=====')
        file = paste0(initial,'_',i,'.rds') # 'combine_','cmv_'
        r = readRDS(paste0(dat.dir,file))
        dim(r) %>% print
        
        res = r %>% select(sample_name,productive_templates) %>% 
            mutate(productive_templates = as.numeric(productive_templates)) %>% 
            dplyr::rename(library.size = productive_templates) %>% unique
        
        saveRDS(res,paste0(output.dir,file))
        
        return(res)
    },mc.cores = detectCores())
    
    concat = do.call(rbind,res.ls)
    saveRDS(concat,paste0(output.dir,initial,'.rds'))
}
