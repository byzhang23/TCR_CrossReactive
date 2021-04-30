## heatmap for Fig2-boxplots
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggpubr)
options(stringsAsFactors = F)

setwd('./covid_adaptive/')
output.dir = './revision/'
font.size = 7

comb = read.csv('./result/manafest_readout/combine_addpseudo.csv')
query = read.csv('./result/manafest_readout/Fig2row.csv') 
query$type = factor(query$type,levels = c("Cross Reactive","COV2S"))
query$sample = factor(query$sample,levels = paste0('CCP',gsub('CCP','',unique(query$sample)) %>% as.numeric %>% sort))
query = query %>% arrange(type,sample,amino_acid)
dat = inner_join(query,comb)
dat$max.ratio = apply(dat[,-c(1:4)],1,function(x) max(x,na.rm = T))
dat$fullname = paste0(dat$amino_acid,'_',dat$sample)
annot = dat[,c('sample','type')]
rownames(annot) = paste0(dat$amino_acid,'_',dat$sample)


sample_colors = c("#F8766D", "#E9842C", "#D69100", "#BC9D00", "#9CA700", "#6FB000", 
                  "#00B813", "#00BD61", "#00C08E", "#00C0B4", "#00BDD4", "#00B5EE", 
                  "#CC79A7", "purple", 
                  "#E69F00","#F0E442", "#D55E00", 
                  "black", "#999999"); 
names(sample_colors) = c('CCP1','CCP2','CCP3','CCP4','CCP5','CCP6','CCP7','CCP8','CCP9','CCP10','CCP11','CCP12',
                         'ED2','ED3',
                         'PC1','PC2','PC3',
                         'JH014_Pre','JH014_Post')
type_colors = c("red", "blue",'darkgreen','pink','yellow'); names(type_colors) = c('Cross Reactive','COV2S','COV2N','NL63S','CEF')
mat_colors = list(sample = sample_colors[names(sample_colors) %in% unique(annot$sample)],
                  type = type_colors[names(type_colors) %in% unique(annot$type)])


## boxplot for crossreactive vs. monoreactive in cov2s
part2 = dat %>% select(fullname,ratio.COV2S,type) %>% rename(group = type,fc = ratio.COV2S)
wilcox.test(part2$fc[part2$group=='Cross Reactive'],part2$fc[part2$group=='COV2S'],alternative = "two")$p.value

png('./revision/Fig3C.png',width = 4,height = 4,res = 300,units = 'in')
ggplot(part2,aes(x=group,y=fc)) +
    geom_boxplot(aes(fill = group),size = 0.8,outlier.size = 0.8) +
    geom_point(size = 0.8) +
    theme_classic() +
    labs(x='',y='Fold Change') +
    stat_compare_means(aes(group = group),label = "p.signif",hide.ns = T,label.y = 160,size = font.size) +
    scale_fill_manual(breaks = c('Cross Reactive','COV2S'),
                      values = c('red','blue')) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 0),
          axis.title = element_text(size = 15), #, face = "bold"
          axis.text = element_text(size = 15), #, face = "bold"
          legend.title=element_text(size=15), 
          legend.text=element_text(size=15),
          plot.title = element_text(size=15)) +
    ylim(0,170)
dev.off()
