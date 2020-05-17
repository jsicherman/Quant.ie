library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(gridExtra)
theme_set(theme_cowplot())

FILE_PATH <- '../output'

if(length(commandArgs(T)) > 0)
  FILE_PATH <- commandArgs(T)[1]

exons <- readRDS(file.path(FILE_PATH, 'quantified', 'counts_exon.rds'))
introns <- readRDS(file.path(FILE_PATH, 'quantified', 'counts_intron.rds'))
qc <- read.csv(file.path(FILE_PATH, 'quantified', 'qc.tsv'), sep = '\t')

# There is a log for every lane (ie. they're unmerged) so we need to do this stupidity.
LOGS <- list.files(file.path(FILE_PATH, 'logs'))

qc.STAR <- do.call(rbind, lapply(names(exons.count), function(sample) {
  do.call('rbind', lapply(LOGS[grepl(sample, LOGS, fixed = T)], function(log) {
    file <- read.delim(file.path(FILE_PATH, 'logs', log), header = F, comment.char = '#', as.is = T)
    
    data.frame(total = file[5, 2] %>% as.integer,
               unique = file[8, 2] %>% as.integer,
               multimapped = file[25, 2] %>% as.integer,
               unmapped = c(file[28, 2], file[30, 2], file[32, 2]) %>% as.integer %>% sum) %>% as.matrix
  })) %>% colSums
})) %>% `rownames<-`(names(exons.count))

exons.count <- Matrix::colSums(exons)
introns.count <- Matrix::colSums(introns)

grid.arrange(data.frame(Sample = names(exons.count),
                        Type = c(rep('Exonic', length(exons.count)), rep('Intronic', length(introns.count)),
                                 rep('Intergenic', length(exons.count)), rep('Multimapped', length(exons.count)),
                                 rep('Unmapped', length(exons.count))),
                        Reads = c(exons.count, introns.count,
                                  qc.STAR[, 'unique'] - exons.count - introns.count, qc.STAR[, 'multimapped'],
                                  qc.STAR[, 'unmapped']) / qc.STAR[, 'total']) %>%
               ggplot(aes(Sample, Reads, fill = Type)) + geom_bar(stat = 'identity') +
               theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = 'none') +
               scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
               scale_x_discrete(breaks = names(exons.count)[seq(1, length(exons.count), by = floor(length(exons.count) / 50) + 1)]) +
               ggtitle('Fractional Counts') + ylab('Counts') + labs(tag = 'A'),
             data.frame(Sample = names(exons.count),
                        Type = c(rep('Exonic', length(exons.count)), rep('Intronic', length(introns.count)),
                                 rep('Intergenic', length(exons.count)), rep('Multimapped', length(exons.count)),
                                 rep('Unmapped', length(exons.count))),
                        Reads = c(exons.count, introns.count,
                                  qc.STAR[, 'unique'] - exons.count, qc.STAR[, 'multimapped'],
                                  qc.STAR[, 'unmapped'])) %>%
               ggplot(aes(Sample, Reads, fill = Type)) + geom_bar(stat = 'identity') +
               theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
               scale_y_continuous(expand = c(0, 0)) +
               scale_x_discrete(breaks = names(exons.count)[seq(1, length(exons.count), by = floor(length(exons.count) / 50) + 1)]) +
               ggtitle('Raw Counts') + ylab(element_blank()) + labs(tag = 'B'),
             nrow = 1) -> grid

ggsave2(file.path(FILE_PATH, 'quantified', 'mapping_distribution.pdf'), grid, width = 11, height = 6)
  
grid.arrange(data.frame(Sample = names(exons.count),
                        Location = c(rep('Exon', length(exons.count)), rep('Intron', length(introns.count))) %>% factor(c('Intron', 'Exon')),
                        Reads = c(exons.count, introns.count) / (exons.count + introns.count)) %>%
               ggplot(aes(Sample, Reads, fill = Location)) + geom_bar(stat = 'identity') +
               theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = 'none') +
               scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
               scale_x_discrete(breaks = names(exons.count)[seq(1, length(exons.count), by = floor(length(exons.count) / 50) + 1)]) +
               ggtitle('Fractional Counts') + ylab('Counts') + labs(tag = 'A'),
             data.frame(Sample = names(exons.count),
                        Location = c(rep('Exon', length(exons.count)), rep('Intron', length(introns.count))) %>% factor(c('Intron', 'Exon')),
                        Reads = c(exons.count, introns.count)) %>%
               ggplot(aes(Sample, Reads, fill = Location)) + geom_bar(stat = 'identity') +
               theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
               scale_y_continuous(expand = c(0, 0)) +
               scale_x_discrete(breaks = names(exons.count)[seq(1, length(exons.count), by = floor(length(exons.count) / 50) + 1)]) +
               ggtitle('Raw Counts') + ylab(element_blank()) + labs(tag = 'B'),
             nrow = 1) -> grid

ggsave2(file.path(FILE_PATH, 'quantified', 'count_distribution.pdf'), grid, width = 11, height = 6)

print('Mitochondrial Reads (%):')
print(setNames(round(qc$mitochondrial / qc$total * 100, 2), rownames(qc)))
print('Ribosomal Reads (%):')
print(setNames(round(qc$rRNA / qc$total * 100, 2), rownames(qc)))