library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(gridExtra)
theme_set(theme_cowplot())

FILE_PATH <- 'output'

if(length(commandArgs(T)) > 0)
  FILE_PATH <- commandArgs(T)[1]

exons <- readRDS(file.path(FILE_PATH, 'exon.rds'))
introns <- readRDS(file.path(FILE_PATH, 'intron.rds'))
qc <- read.csv(file.path(FILE_PATH, 'qc.tsv'), sep = '\t') %>% `rownames<-`(.[, 1]) %>% .[, -1]

exons.count <- Matrix::colSums(exons)
introns.count <- Matrix::colSums(introns)

grid.arrange(data.frame(Sample = names(exons.count),
                        Location = c(rep('Exon', length(exons.count)), rep('Intron', length(introns.count))) %>% factor(c('Intron', 'Exon')),
                        Reads = c(exons.count, introns.count) / (exons.count + introns.count)) %>%
               ggplot(aes(Sample, Reads, fill = Location)) + geom_bar(stat = 'identity') +
               theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
               scale_y_continuous(expand = c(0, 0)) + ggtitle('Fractional Counts') + ylab('Counts') + labs(tag = 'A'),
             data.frame(Sample = names(exons.count),
                        Location = c(rep('Exon', length(exons.count)), rep('Intron', length(introns.count))) %>% factor(c('Intron', 'Exon')),
                        Reads = c(exons.count, introns.count)) %>%
               ggplot(aes(Sample, Reads, fill = Location)) + geom_bar(stat = 'identity') +
               theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
               scale_y_continuous(expand = c(0, 0)) + ggtitle('Raw Counts') + ylab(element_blank()) + labs(tag = 'B'),
             nrow = 1) -> grid

ggsave2(file.path(FILE_PATH, 'count_distribution.pdf'), grid, width = 11, height = 6)

print('Mitochondrial Reads (%):')
print(setNames(round(qc$mitochondrial / qc$total * 100, 2), rownames(qc)))
print('Ribosomal Reads (%):')
print(setNames(round(qc$rRNA / qc$total * 100, 2), rownames(qc)))