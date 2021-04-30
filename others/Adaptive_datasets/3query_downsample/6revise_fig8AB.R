## barplots (Fig8, cross reactive clones exist in either query datasets)
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


bar<- function(target = '1CrossReactive',
                         ht = 15,
                         barplot = T,
                         font.size = 7){

    
    ## percent of patients
    ttl.covid = 1477
    ttl.cmv = 709
    comb1.ls = readRDS(paste0(downsample.dir,target,'_summary.rds'))
    draw = do.call(rbind,lapply(1:length(comb1.ls),function(i){
        d = comb1.ls[[i]]
        d$cmv = d$cmv/ttl.cmv
        d$covid = d$covid/ttl.covid
        keep.id = which(rowSums(d[,2:3],na.rm = T) > 0)
        
        tmp = d[keep.id,] %>% mutate(sample_name = names(comb1.ls)[i]) %>% 
         #   filter(amino_acid %in% dat$amino_acid) %>% 
            gather(type,percent,-c(amino_acid,sample_name)) %>%
            mutate(type = ifelse(type=='covid','Covid','Pre-Covid'))
        
        tmp[is.na(tmp)] = 0
        tmp
    }))
    draw$type = factor(draw$type, levels = c('Covid','Pre-Covid'))
    draw$amino_acid = factor(draw$amino_acid,levels = unique(draw$amino_acid) %>% sort)
    
    png(paste0('./revision/Fig8A_barplot_',target,'_fdr.png'),width = 6,height = ht,units = 'in',res = 300) # two plots: 14
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
        #grid.arrange(p1,p2,ncol=2,widths=c(6,7))
        dev.off()
        
    ## Fig8B boxplots
        png(paste0('./revision/Fig8B_boxplot_',target,'.png'),width = 6,height = 6,units = 'in',res = 300) # two plots: 14
        set.seed(12345)
        b <- runif(nrow(draw), -0.02, 0.02)
        draw$type = factor(draw$type,levels = c('Pre-Covid','Covid'))
        
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
