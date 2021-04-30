## heatmap/dendrogram (for query ones)
library(dplyr)
library(stringdist) 
library(data.table)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
options(stringsAsFactors = F)
setwd('./covid_adaptive/')

block_size <- function(new,
                       min.dist){
    
    block = c()
    for(bs in c(2:10)){ # general usage: replace 10 by nrow(new)
        
        max.value = c()
        last = nrow(new)-bs + 1
        for(i in 1:last){
            
            max.value = c(max.value,max(new[i:(i+bs-1),i:(i+bs-1)]))
            
        }
        lowest = min(max.value) 
        block = c(block,ifelse(lowest<=min.dist,1,0))
    }
    
    return(sum(block)+1)
    
}

findblock <- function(files = NULL,
                 input = NULL,
                 full = F,
                 ht = 10,
                 min.dist = 3,
                 plotout = F){
    
    
    if(!is.null(files)){
        dat = do.call(rbind,lapply(files,function(f){
            
            read.csv(paste0('./revision/',f,'.csv')) 
        }))
        dat$type = sub('^[0-9]','',dat$type)
    }else{
        dat = input
        dat$type = sub('^[0-9]','',dat$type)
    }
    
    
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
    type_colors = c("red", "blue",'darkgreen','pink','lightgreen','lightblue','grey'); names(type_colors) = c('CrossReactive','COV2S_Specific','COV2N_Specific','NL63S_Specific','CEF','EBV','HIV')
    mat_colors = list(Pt.ID = sample_colors[names(sample_colors) %in% unique(annotation_col$Pt.ID)],
                      'TCR reactivity' = type_colors[names(type_colors) %in% unique(annotation_col$`TCR reactivity`)])
    
    if(plotout==T){
        ## pheatmap
        pdf(paste0('./revision/',files,'.pdf'),height = ht,width = ht)
        # breaksList = c(0:14,seq(16, max(mat), by = 10))
        breaksList = c(0,1,2,3,max(mat))
        ht = pheatmap(mat,
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
    }else{
        breaksList = c(0,1,2,3,max(mat))
        ht = pheatmap(mat,
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
    }
    
    
    rows = ht$tree_row$labels[ht$tree_row$order]
    new = mat[rows,rows]
    
    
    
    res = block_size(new,min.dist)
}


## example
findblock(files = c('HIV_seed1'),ht = 15,full = F,min.dist = 3,plotout = T)

## calculate block size
cr = read.csv('./data/query_revision/1CrossReactive.csv')
info = cr %>% group_by(sample) %>% summarise(n = n_distinct(amino_acid))
dat.ls = readRDS('./revision/HIV_A.rds')
sampling <- function(dat.ls,
                     info,
                     seed = 1){
    
    
    samp.ls = lapply(dat.ls,function(dt){
        
        n = info$n[info$sample==dt$sample %>% unique]
        set.seed(seed)
        sub = dt[sample(1:nrow(dt),n,replace = F),]
        
        sub %>% select(-freq) %>% rename(amino_acid = aminoAcid) %>% mutate(type = 'HIV') %>% select(type,sample,amino_acid)
    })
    
    do.call(rbind,samp.ls)
}

bl = c()
for(s in 1:100){
    
    dat.random = sampling(dat.ls,info,seed = s)
    res = findblock(files = NULL,input = dat.random,ht = 15,full = F,min.dist = 3,plotout = F)
    bl = c(bl,res)
}
write.csv(bl,'./revision/HIV_blocksize.csv',row.names = F)

bl = data.frame(block_size = bl)

pdf('./revision/HIV_blocksize.pdf',height = 3,width = 5)
ggplot(bl,aes(x=block_size)) +
    geom_histogram() +
    xlim(0,10) +
    labs(y='Counts') +
    geom_vline(xintercept = 10,color = 'red',linetype = 'dashed') +
    theme_bw() +
    scale_x_continuous(breaks=c(0:10))
dev.off()
