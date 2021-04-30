## heatmap zoom in
library(dplyr)
library(stringdist) 
library(data.table)
library(pheatmap)
library(RColorBrewer)
options(stringsAsFactors = F)
setwd('./covid_adaptive/')


heat_zoom <- function(files = c('2COV2S_Specific'),
                      full = F,
                      from = '', # data
                      to = '', # data
                      ht = 5){
    
    
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
    
    breaksList = seq(min(mat), max(mat), by = 1)
    pt = pheatmap(mat,
                  color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(length(breaksList)),
                  breaks = breaksList,
                  fontsize = 8,
                  cellwidth = 8,cellheight = 8,
                  border_color = NA,
                  annotation_colors = mat_colors,
                  annotation_col = annotation_col,
                  annotation_row = annotation_col,
                  show_rownames = T,
                  show_colnames = F)
    
    order.label = colnames(mat[,pt$tree_col[["order"]]])
    
    ord = mat[order.label,order.label]
    idx = seq(grep(from,order.label),grep(to,order.label),1)
    sub = ord[idx,idx]
    
    png(paste0('./revision/heatmap-zoomin_',paste(files,collapse = '_'),from,'.png'),height = ht,width = ht,units = 'in',res = 300)
    pheatmap(sub,
             color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(length(breaksList)),
             breaks = breaksList,
             fontsize = 8,
             cellwidth = 8,cellheight = 8,
             border_color = NA,
             annotation_colors = mat_colors,
             annotation_col = annotation_col[rownames(annotation_col) %in% rownames(sub),],
             annotation_row = annotation_col[rownames(annotation_col) %in% rownames(sub),],
             cluster_rows = F,
             cluster_cols = F,
             angle_col = '90',
             show_rownames = T,
             show_colnames = F)
    dev.off()
}
