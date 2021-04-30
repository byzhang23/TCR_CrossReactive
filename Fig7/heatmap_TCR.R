## heatmap/dendrogram (for query ones)
library(dplyr)
library(stringdist) 
library(data.table)
library(pheatmap)
library(RColorBrewer)
options(stringsAsFactors = F)
setwd('./covid_adaptive/')

files = c('1CrossReactive','2COV2S_Specific')
heat <- function(files = c('1CrossReactive','2COV2S_Specific'),
                 full = F,
                 ht = 10){
    
    dat = do.call(rbind,lapply(files,function(f){
        
        read.csv(paste0('./data/query_revision/',f,'.csv')) 
    }))
    dat$type = sub('^[0-9]','',dat$type)
    
    ## remove first 3 and last 3 AAs
    dat$middle = substr(dat$amino_acid,start = 4,stop = nchar(dat$amino_acid)-3)
    
    ## distance
    if(full){
        mat = do.call(rbind,lapply(dat$amino_acid,function(q){
            
            stringdist(q,dat$amino_acid,method = 'lv')
        }))
        rownames(mat) = colnames(mat) = paste0(dat$amino_acid,'_',dat$sample)
    }else{
        mat = do.call(rbind,lapply(dat$middle,function(q){
            
            stringdist(q,dat$middle,method = 'lv')
        }))
        rownames(mat) = colnames(mat) = paste0(dat$amino_acid,'_',dat$sample)
        
    }
    
    ## annotation
    annotation_col = data.frame(Pt.ID = dat$sample,
                                'TCR reactivity' = dat$type,check.names = F)
    rownames(annotation_col) = rownames(mat)
    sample_colors = c("#F8766D", "#E9842C", "#D69100", "#BC9D00", "#9CA700", "#6FB000", 
                      "#00B813", "#00BD61", "#00C08E", "#00C0B4", "#00BDD4", "#00B5EE", 
                      "#7F96FF", 'orchid',"purple", "#E69F00","#F0E442", "#D55E00","firebrick", 
                      "black", "#999999"); 
    names(sample_colors) = c('CCP1','CCP2','CCP3','CCP4','CCP5','CCP6','CCP7','CCP8','CCP9','CCP10','CCP11','CCP12',
                             'PC1','PC2','PC3','PC4','PC5','PC6','PC7',
                             'ED2','ED3')
    type_colors = c("red", "blue",'darkgreen','pink','lightgreen','lightblue'); names(type_colors) = c('CrossReactive','COV2S_Specific','COV2N_Specific','NL63S_Specific','CEF','EBV')
    mat_colors = list(Pt.ID = sample_colors[names(sample_colors) %in% unique(annotation_col$Pt.ID)],
                      'TCR reactivity' = type_colors[names(type_colors) %in% unique(annotation_col$`TCR reactivity`)])
    
    ## pheatmap
    pdf(paste0('./revision/heatmap_',paste(files,collapse = '_'),'.pdf'),height = ht,width = ht)
    breaksList = seq(min(mat), max(mat), by = 1)
    pheatmap(mat,
                  color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(length(breaksList)),
                  breaks = breaksList,
                  fontsize = 8,
                  cellwidth = 8,cellheight = 8,
                  border_color = NA,
                  annotation_colors = mat_colors,
                  annotation_col = annotation_col,
                  annotation_row = annotation_col,
                    angle_col = '90',
                  show_rownames = T,
                  show_colnames = F)
    dev.off()
    
    png(paste0('./revision/heatmap_',paste(files,collapse = '_'),'.png'),height = ht,width = ht,units = 'in',res = 300)
    breaksList = seq(min(mat), max(mat), by = 1)
    pheatmap(mat,
             color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(length(breaksList)),
             breaks = breaksList,
             fontsize = 8,
             cellwidth = 8,cellheight = 8,
             border_color = NA,
             annotation_colors = mat_colors,
             annotation_col = annotation_col,
             annotation_row = annotation_col,
             angle_col = '90',
             show_rownames = T,
             show_colnames = F)
    dev.off()
    
}

## example
heat(files = c('1CrossReactive'),ht = 12,full = F)
heat(files = c('2COV2S_Specific'),ht = 34,full = F)
heat(files = c('3COV2N_Specific'),ht = 32,full = F)
heat(files = c('4CEF'),ht = 44,full = F)
