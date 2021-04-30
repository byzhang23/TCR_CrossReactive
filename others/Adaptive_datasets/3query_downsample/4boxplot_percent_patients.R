## boxplot/volin plots
library(ggplot2)
library(data.table)
library(ggpubr)
library(tidyr)
library(dplyr)
library(gridExtra)
options(stringsAsFactors = F)

setwd('./covid_adpative/')
cmv.dir = './result/query_cmv/'
covid.dir = './result/query/'
downsample.dir = './revision/'


target = '1CrossReactive'

compare_patients <- function(target = '1CrossReactive',
                wd = 8,
                ht = 5){
    
    dat1 = fread(paste0(covid.dir,target,'_detail.csv')) %>% select(sample_name,amino_acid,templates,productive_frequency) %>% mutate(type = 'Covid')
    dat2 = fread(paste0(cmv.dir,target,'_detail_cmv.csv')) %>% mutate(type = 'Pre-Covid')
    comb = rbind(dat1,dat2)
    common = intersect(unique(dat1$amino_acid),unique(dat2$amino_acid))
    
    dat = comb %>% filter(amino_acid %in% common)
    dat$type = factor(dat$type, levels = c('Covid','Pre-Covid'))
    dat$amino_acid = factor(dat$amino_acid,levels = unique(dat$amino_acid) %>% sort)
    
    
    ## percent of patients
    ttl.covid = 1477
    ttl.cmv = 709
    comb1.ls = readRDS(paste0(downsample.dir,target,'_summary.rds'))
    draw = do.call(rbind,lapply(1:length(comb1.ls),function(i){
        d = comb1.ls[[i]]
        d$cmv = d$cmv/ttl.cmv
        d$covid = d$covid/ttl.covid
        tmp = d %>% mutate(sample_name = names(comb1.ls)[i]) %>% 
            filter(amino_acid %in% dat$amino_acid) %>% 
            gather(type,percent,-c(amino_acid,sample_name)) %>%
            mutate(type = ifelse(type=='covid','Covid','Pre-Covid'))
        
        tmp[is.na(tmp)] = 0
        tmp
    }))
    
    # median
    draw = draw %>% group_by(amino_acid,type) %>% summarise(percent = median(percent))
    draw$type = factor(draw$type, levels = c('Covid','Pre-Covid'))
    draw$amino_acid = factor(draw$amino_acid,levels = unique(draw$amino_acid) %>% sort)
    print(identical(levels(draw$amino_acid),levels(dat$amino_acid)))
    
    
        set.seed(12345)
        b <- runif(nrow(draw), -0.02, 0.02)
        draw$type = factor(draw$type,levels = c('Pre-Covid','Covid'))
        
        png(paste0('./revision/boxplot_percentpatients_',target,'.png'),width = wd,height = ht,units = 'in',res = 300)
        (ggplot(draw,aes(x = as.numeric(type),y=percent,fill = type)) +
                geom_boxplot(outlier.shape = NA,width = 0.7) +
                geom_point(aes(x = as.numeric(type) + b, y = percent),size = 0.6) +
                geom_line(aes(x  = as.numeric(type) + b, y = percent, group = amino_acid),color = 'grey',alpha = 0.5) +
                #geom_line(aes(group = amino_acid),color = 'grey') +
                stat_compare_means(aes(group = type),paired = T, label = "p.signif",hide.ns = T) + 
                
                theme_classic() +
                scale_x_continuous(breaks = c(1,2), labels = c("Pre-Covid","Covid"))+
                scale_fill_manual(breaks = c('Pre-Covid','Covid'),
                                  values = c('orange','red')) +
                labs(x='',y='Percent of Patients (Median)') +
                scale_y_continuous(labels=scales::percent) +
                theme(axis.text.x = element_text(angle = 0),
                      legend.position = 'none',
                      axis.title = element_text(size = 15), #, face = "bold"
                      axis.text = element_text(size = 15), #, face = "bold"
                      legend.title=element_text(size=15), 
                      legend.text=element_text(size=15),
                      plot.title = element_text(size=15)) ) %>% print
        dev.off()
    
}

## example
## 1. Barplot
compare_patients(target = '1CrossReactive',wd = 8,ht = 8,bar = T)
bar(target = '2COV2S_Specific',wd = 12,ht = 20,bar = T)
bar(target = '3COV2N_Specific',wd = 8,ht = 12,bar = T)
#bar(target = '4NL63S_Specific',wd = 8,ht = 12)


## 2. Boxplot
compare_patients(target = '1CrossReactive',wd = 4,ht = 4)
compare_patients(target = '2COV2S_Specific',wd = 4,ht = 4)
compare_patients(target = '3COV2N_Specific',wd = 4,ht = 4)
