## sequence logo
suppressMessages(library(motifStack))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
library(Biostrings)
options(stringsAsFactors = F)


output.dir = ''
dat = read.csv(paste0(output.dir,'sequence_logo.csv')) %>% mutate(length = nchar(cross_reactive))
mat = AAStringSet(dat$cross_reactive) %>% consensusMatrix(as.prob=T) %>% as.matrix


png(paste0(output.dir,"sequence_logo.png"),res = 300,width = 8,height = 5,units = 'in')
pfm1 = new("pfm", mat=mat,name = '',color=colorset(alphabet="AA",colorScheme="chemistry"))
plot(pfm1, ic.scale=FALSE, ylab="probability")
dev.off() 
