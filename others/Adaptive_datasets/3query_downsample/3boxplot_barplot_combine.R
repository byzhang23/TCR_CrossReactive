## boxplot/volin plots
library(ggplot2)
library(data.table)
library(ggpubr)
library(dplyr)
library(tidyr)
library(gridExtra)
options(stringsAsFactors = F)

setwd('./covid_adaptive/')
cmv.dir = './result/query_cmv/'
covid.dir = './result/query/'
downsample.dir = './revision/'

target = '1CrossReactive'

box_bar_comb <- function(target = '1CrossReactive',
                         ht = 15,
                         barplot = T,
                         font.size = 7){
    
    ## frequency
    dat1 = fread(paste0(covid.dir,target,'_detail.csv')) %>% select(sample_name,amino_acid,templates,productive_frequency) %>% mutate(type = 'Covid')
    dat2 = fread(paste0(cmv.dir,target,'_detail_cmv.csv')) %>% mutate(type = 'Pre-Covid')
    comb = rbind(dat1,dat2)
    common = intersect(unique(dat1$amino_acid),unique(dat2$amino_acid))
    
    dat = comb %>% filter(amino_acid %in% common)
    dat$type = factor(dat$type, levels = c('Covid','Pre-Covid'))
    dat$amino_acid = factor(dat$amino_acid,levels = unique(dat$amino_acid) %>% sort)
    
    comp.res = compare_means(productive_frequency~type,data = dat,group.by = 'amino_acid',p.adjust.method = 'BH')
    write.csv(comp.res,paste0('./revision/pvalue_frequency_',target,'.csv'),row.names = F)
    
    
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
    draw$type = factor(draw$type, levels = c('Covid','Pre-Covid'))
    draw$amino_acid = factor(draw$amino_acid,levels = unique(draw$amino_acid) %>% sort)
    print(identical(levels(draw$amino_acid),levels(dat$amino_acid)))
    
    pp.res = compare_means(percent~type,data = draw,group.by = 'amino_acid',p.adjust.method = 'BH')
    write.csv(pp.res,paste0('./revision/pvalue_percentpatients_',target,'.csv'),row.names = F)
    
    if(!barplot){
        png(paste0('./revision/combine_',target,'_fdr.png'),width = 14,height = ht,units = 'in',res = 300) # two plots: 14
        
        # change p.signif
        comp.res$p.signif = ifelse(comp.res$p.adj>0.05,'ns',
                                   ifelse(comp.res$p.adj <= 0.0001,'****',
                                          ifelse(comp.res$p.adj <=  0.001,'***',
                                                 ifelse(comp.res$p.adj <= 0.01,'**','*'))))
        
        pp.res$p.signif = ifelse(pp.res$p.adj>0.05,'ns',
                                   ifelse(pp.res$p.adj <= 0.0001,'****',
                                          ifelse(pp.res$p.adj <=  0.001,'***',
                                                 ifelse(pp.res$p.adj <= 0.01,'**','*'))))
        
        p2 = (ggplot(dat,aes(x=amino_acid,y=log10(productive_frequency*100))) +
                  geom_boxplot(aes(fill = type),outlier.size = 0.6,width = 0.7) +
                  stat_pvalue_manual(comp.res, label = "p.signif", x='amino_acid',y.position = 0,hide.ns = T,size = font.size,color = "red") +
                  
                  theme_classic() +
                  scale_fill_manual(breaks = c('Pre-Covid','Covid'),
                                    values = c('orange','red')) +
                  labs(x='',y='log10(Productive Frequency*100%)') +
                  coord_flip() +
                  #theme(axis.text.y = element_blank()) +
                  theme(axis.text.x = element_text(angle = 0),
                        axis.title = element_text(size = 15), #, face = "bold"
                        axis.text = element_text(size = 15), #, face = "bold"
                        legend.title=element_text(size=15), 
                        legend.text=element_text(size=15),
                        plot.title = element_text(size=15),
                        legend.position = 'right') +
                  theme(axis.text.y = element_blank(), 
                        #axis.ticks.y = element_blank(), 
                        axis.title.y = element_blank())) %>% print
        p2$layers[[2]]$aes_params$textsize <- font.size
        
        
        p1 = (ggplot(draw,aes(x=amino_acid,y=percent)) +
                  geom_boxplot(aes(fill = type),outlier.size = 0.6,width = 0.7) +
                  # stat_pvalue_manual(pp.res, label = "p.signif", x='amino_acid',y.position = max(draw$percent) + 0.1,hide.ns = T,size = font.size,color = "red") +
                  
                  theme_classic() +
                  scale_fill_manual(breaks = c('Pre-Covid','Covid'),
                                    values = c('orange','red')) +
                  labs(x='',y='Percent of Patients') +
                  scale_y_continuous(labels=scales::percent) +
                  coord_flip() +
                  theme(axis.text.x = element_text(angle = 0),
                        axis.title = element_text(size = 15), #, face = "bold"
                        axis.text = element_text(size = 15), #, face = "bold"
                        legend.title=element_text(size=15), 
                        legend.text=element_text(size=15),
                        plot.title = element_text(size=15),
                        legend.position = 'none') ) %>% print
        # p1$layers[[2]]$aes_params$textsize <- font.size
        
        # barplot
        # p1 = (ggplot(draw,aes(x=amino_acid,y=percent,fill = type)) +
        #           geom_bar(stat='identity',width=.7, position = "dodge") +
        #           #geom_violin() +
        #           stat_compare_means(aes(group = type), label = "p.signif",hide.ns = T) + 
        #           
        #           theme_classic() +
        #           scale_fill_manual(breaks = c('Pre-Covid','Covid'),
        #                             values = c('orange','red')) +
        #           labs(x='',y='Percent of Patients') +
        #           scale_y_continuous(labels=scales::percent) +
        #           coord_flip() +
        #           theme(axis.text.x = element_text(angle = 0),
        #                 legend.position = 'none',
        #                 axis.title = element_text(size = 15), #, face = "bold"
        #                 axis.text = element_text(size = 15), #, face = "bold"
        #                 legend.title=element_text(size=15), 
        #                 legend.text=element_text(size=15),
        #                 plot.title = element_text(size=15)) ) %>% print
        grid.arrange(p1,p2,ncol=2,widths=c(6,7))
        dev.off()
        
    }else{
        png(paste0('./revision/combine_bar_',target,'_fdr.png'),width = 14,height = ht,units = 'in',res = 300) # two plots: 14
        
        # change p.signif
        comp.res$p.signif = ifelse(comp.res$p.adj>0.05,'ns',
                                   ifelse(comp.res$p.adj <= 0.0001,'****',
                                          ifelse(comp.res$p.adj <=  0.001,'***',
                                                 ifelse(comp.res$p.adj <= 0.01,'**','*'))))
        
        pp.res$p.signif = ifelse(pp.res$p.adj>0.05,'ns',
                                 ifelse(pp.res$p.adj <= 0.0001,'****',
                                        ifelse(pp.res$p.adj <=  0.001,'***',
                                               ifelse(pp.res$p.adj <= 0.01,'**','*'))))
        
        p2 = (ggplot(dat,aes(x=amino_acid,y=log10(productive_frequency*100))) +
                  geom_boxplot(aes(fill = type),outlier.size = 0.6,width = 0.7) +
                  stat_pvalue_manual(comp.res, label = "p.signif", x='amino_acid',y.position = 0,hide.ns = T,size = font.size,color = "red") +
                  
                  theme_classic() +
                  scale_fill_manual(breaks = c('Pre-Covid','Covid'),
                                    values = c('orange','red')) +
                  labs(x='',y='log10(Productive Frequency*100%)') +
                  coord_flip() +
                  #theme(axis.text.y = element_blank()) +
                  theme(axis.text.x = element_text(angle = 0),
                        axis.title = element_text(size = 15), #, face = "bold"
                        axis.text = element_text(size = 15), #, face = "bold"
                        legend.title=element_text(size=15), 
                        legend.text=element_text(size=15),
                        plot.title = element_text(size=15),
                        legend.position = 'right') +
                  theme(axis.text.y = element_blank(), 
                        #axis.ticks.y = element_blank(), 
                        axis.title.y = element_blank())) %>% print
        p2$layers[[2]]$aes_params$textsize <- font.size
        
        draw = draw %>% group_by(amino_acid,type) %>% summarise(percent = median(percent))
        # barplot
        p1 = (ggplot(draw,aes(x=amino_acid,y=percent,fill = type)) +
                  geom_bar(stat='identity',width=.7, position = "dodge") +
                  #geom_violin() +
                  stat_compare_means(aes(group = type), label = "p.signif",hide.ns = T) +

                  theme_classic() +
                  scale_fill_manual(breaks = c('Pre-Covid','Covid'),
                                    values = c('orange','red')) +
                  labs(x='',y='Percent of Patients (Median)') +
                  scale_y_continuous(labels=scales::percent) +
                  coord_flip() +
                  theme(axis.text.x = element_text(angle = 0),
                        legend.position = 'none',
                        axis.title = element_text(size = 15), #, face = "bold"
                        axis.text = element_text(size = 15), #, face = "bold"
                        legend.title=element_text(size=15),
                        legend.text=element_text(size=15),
                        plot.title = element_text(size=15)) ) %>% print
        grid.arrange(p1,p2,ncol=2,widths=c(6,7))
        dev.off()
        
        
    }
    
}

## example

## barplot
box_bar_comb(target = '1CrossReactive', ht = 10, barplot = T)
box_bar_comb(target = '2COV2S_Specific', ht = 30, barplot = T)
box_bar_comb(target = '3COV2N_Specific', ht = 25, barplot = T)

## boxplots (most of them are identical across 100 downsampling experiment)
box_bar_comb(target = '1CrossReactive', ht = 10, barplot = F)
box_bar_comb(target = '2COV2S_Specific', ht = 30, barplot = F)
box_bar_comb(target = '3COV2N_Specific', ht = 25, barplot = F)

