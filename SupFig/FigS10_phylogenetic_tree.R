## heatmap/dendrogram (for query ones)
library(dplyr)
library(stringdist) 
library(data.table)
library(pheatmap)
library(ape) #polygenetic tree
options(stringsAsFactors = F)
setwd('./covid_adaptive/')

files = c('1CrossReactive','2COV2S_Specific')
dendro <- function(files = c('1CrossReactive','2COV2S_Specific'),
                 full = F,
                 ht = 10,
                 cex = 0.2){
    
    dat = do.call(rbind,lapply(files,function(f){
        
        read.csv(paste0('./data/query_revision/',f,'.csv')) 
    }))
    dat$type = sub('^[0-9]','',dat$type)
    dat$nam = paste0(dat$amino_acid,'...',dat$sample)
    
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
        rownames(mat) = colnames(mat) = paste0(dat$amino_acid,'...',dat$sample)
        
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
    mat_colors = data.frame(sample.col = sample_colors[names(sample_colors) %in% unique(annotation_col$Pt.ID)])
    mat_colors$sample = rownames(mat_colors)
    
    
    ## polytree
    phylotree <- function(res = res,
                          type = 'unrooted'){
        
        hc = hclust(dist(res))
        poly.tree = as.phylo(hc)
        
        if(!is.null(dat)){
            
            annot = left_join(left_join(data.frame(nam = poly.tree$tip.label),dat),mat_colors)
            plot(poly.tree, type = type, cex = cex,no.margin = F,font = 2, tip.color = annot$sample.col,label.offset = 1,lab4ut = 'axial')
        }else{
            
            plot(poly.tree, type = type, cex = cex,no.margin = F,font = 2)
        }
        
        
        
    }
    
    
    ## pheatmap
    png(paste0('./revision/phylo_',paste(files,collapse = '_'),'.png'),height = ht,width = ht,units = 'in',res = 300)
    phylotree(mat,'unrooted')
    dev.off()
    
}

## example
dendro(files = c('1CrossReactive'),ht = 15,full = F,cex = 0.8)
