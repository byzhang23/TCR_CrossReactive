## heatmap for Fig2
library(dplyr)
library(data.table)
library(RColorBrewer)
library(tidyr)
library(pheatmap)
options(stringsAsFactors = F)
setwd('./covid_adaptive/')

output.dir = './revision/'

comb = read.csv('./result/manafest_readout/combine_addpseudo.csv')
query = read.csv('./result/manafest_readout/Fig2row.csv') 
query$type = factor(query$type,levels = c("Cross Reactive","COV2S"))
query$sample = factor(query$sample,levels = paste0('CCP',gsub('CCP','',unique(query$sample)) %>% as.numeric %>% sort))
query = query %>% arrange(type,sample,amino_acid)
dat = inner_join(query,comb)

mat = dat[,grep('^ratio',colnames(dat))]
rownames(mat) = paste0(dat$amino_acid,'_',dat$sample)
colnames(mat) = gsub('ratio.','',colnames(mat))
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

## order
# sub1.id = rownames(annot)[annot$type=='Cross Reactive']
# sub2.id = rownames(annot)[annot$type=='COV2S']
# sub1 = mat[rownames(mat) %in% sub1.id,]
# sub1[is.na(sub1)] = 0
# sub2 = mat[rownames(mat) %in% sub2.id,]
# sub2[is.na(sub2)] = 0
# hc1 = hclust(dist(sub1))
# label1 = hc1$labels[hc1$order]
# hc2 = hclust(dist(sub2))
# label2 = hc2$labels[hc2$order]
# labels = c(label1,label2)


# breaksList = c(seq(0, 100, by = 5),seq(100,max(mat,na.rm=T),50)) %>% unique
breaksList = c(seq(0, 100, by = 5),max(mat,na.rm=T)) %>% unique
# pdf(paste0(output.dir,'Fig3A_main.pdf'),width = 9,height = 10)
# pheatmap(#mat[labels,],
#     mat[order(mat[,1],decreasing = T),], # order by descending fold change
#     color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaksList)),
#     breaks = breaksList,
#     fontsize = 7,
#     #cellwidth = 8,cellheight = 8,
#     border_color = 'black',
#     annotation_colors = mat_colors,
#     annotation_row = annot,
#     cluster_cols = F,
#     cluster_rows = F,
#     angle_col = 45,
#     show_rownames = T,
#     show_colnames = T)
# dev.off()

## png
png(paste0(output.dir,'Fig3A_main.png'),width = 9,height = 10,res = 300,units = 'in')
ht = pheatmap(#mat[labels,],
    mat[order(mat[,1],decreasing = T),], # order by descending fold change
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaksList)),
    breaks = breaksList,
    fontsize = 7,
    #cellwidth = 8,cellheight = 8,
    border_color = 'black',
    annotation_colors = mat_colors,
    annotation_row = annot,
    cluster_cols = F,
    cluster_rows = F,
    angle_col = 45,
    show_rownames = T,
    show_colnames = T)
dev.off()

aa = mat[order(mat[,1],decreasing = T),]
saveRDS(rownames(aa),'./revision/Fig3A_main_roworder.rds')

