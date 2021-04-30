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


## boxplot for COV2S vs. others in Crossreactive (old)
part1 = dat %>% filter(type=='Cross Reactive') %>%
    select(fullname,ratio.COV2S,max.ratio) %>% rename(COV2S = ratio.COV2S,Others= max.ratio) %>% 
    gather(group,fc,-c(fullname))

wilcox.test(part1$fc[part1$group=='COV2S'],part1$fc[part1$group=='Others'],alternative = "two",paired=T)$p.value

png('./revision/Fig3B_supplemental.png',width = 4,height = 4,res = 300,units = 'in')
ggplot(part1,aes(x=group,y=fc)) +
    geom_boxplot(aes(fill = group),size = 0.8,outlier.size = 0.8) +
    geom_point(size = 0.8) +
    geom_line(aes(group = fullname),color = 'grey',alpha = 0.6) +
    theme_classic() +
    labs(x='',y='Fold Change') +
    stat_compare_means(aes(group = group),paired = T,label = "p.signif",hide.ns = T,label.y = 380,size = font.size) +
    scale_fill_manual(breaks = c('COV2S','Others'),
                      values = c('purple','yellow')) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 0),
          axis.title = element_text(size = 15), #, face = "bold"
          axis.text = element_text(size = 15), #, face = "bold"
          legend.title=element_text(size=15), 
          legend.text=element_text(size=15),
          plot.title = element_text(size=15))
dev.off()

# ## 1st try: shuffle among viruses and repeat
# main = dat[which(dat$type=='Cross Reactive'),paste0('ratio.',c('COV2S','NL63S','229ES','HKU1S','OC43S'))]
# library(combinat)
# permute.all = do.call(rbind,permn(1:ncol(main)))[-1,]
# 
# p.ls = c()
# len = c()
# mean.ls = c()
# perm.mat = c()
# for(s in 1:nrow(permute.all)){
#     
#     shuffle = permute.all[s,]
#     new = main[,shuffle]
#     colnames(new) = colnames(main)
#     ratio.COV2S = new[,1]
#     max.ratio = apply(new[,-1],1,function(x) max(x,na.rm = T))
#     diff = max.ratio - ratio.COV2S
#     perm.mat = cbind(perm.mat,diff)
#     
#     
#     mean.ls = c(mean.ls,mean(diff,na.rm = T))
#     len = c(len,sum(!is.na(diff)))
#     p.ls = c(p.ls,wilcox.test(diff,alternative = 'greater')$p.value)
# 
# }
# 
# obs = dat$max.ratio[dat$type=='Cross Reactive'] - dat$ratio.COV2S[dat$type=='Cross Reactive']
# comb.mat = cbind(obs,perm.mat)
# colnames(comb.mat)=c('obs',paste0('permute_',1:nrow(permute.all)))
# pvalue = apply(comb.mat,1,function(x) sum(x[-1]>x[1],na.rm=T))/nrow(permute.all)
# fdr = p.adjust(pvalue,method = 'BH')
# comb.dat = comb.mat %>% data.frame %>% mutate(pvalue = pvalue,fdr = fdr,fullname = dat$fullname[dat$type=='Cross Reactive'])
# write.csv(comb.dat[,c(123,1,121:122,2:120)],paste0(output.dir,'Fig3B_permute.csv'),row.names = F)
# 




## 2nd try: pairwise comparison (still paired one sided test, but sample size differs across comparisons)
# paired one-sided test (hypothesize for cross reactive)
pairwise = dat %>% filter(type=='Cross Reactive') %>% 
    mutate(NL63S = ratio.NL63S - ratio.COV2S,
           '229ES' = ratio.229ES - ratio.COV2S,
           HKU1S = ratio.HKU1S - ratio.COV2S,
           OC43S = ratio.OC43S - ratio.COV2S)

pair.summary = do.call(rbind,
                        lapply(c('NL63S','229ES','HKU1S','OC43S'),function(x){
                            
                            sub = pairwise[!is.na(pairwise[,x]),]
                            sample.size = nrow(sub)
                            p = wilcox.test(sub[,x],alternative = 'greater')$p.value #others are higher than cov2s
                            #p = wilcox.test(sub[,x],alternative = 'two')$p.value # difference
                            data.frame(peptide = x, sample.size = sample.size,pvalue = p)
                        }))



