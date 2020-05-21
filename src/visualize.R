library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(reshape2)
theme_set(theme_cowplot())

FILE_PATH <- '../output'

if(length(commandArgs(T)) > 0)
  FILE_PATH <- commandArgs(T)[1]

exons <- readRDS(file.path(FILE_PATH, 'quantified', 'counts_exon.rds'))
introns <- readRDS(file.path(FILE_PATH, 'quantified', 'counts_intron.rds'))
qc <- read.csv(file.path(FILE_PATH, 'quantified', 'qc.tsv'), sep = '\t')

exons.count <- Matrix::colSums(exons)
introns.count <- Matrix::colSums(introns)

# There is a log for every lane (ie. they're unmerged) so we need to do this stupidity.
LOGS <- list.files(file.path(FILE_PATH, 'logs'))

qc.STAR <- do.call(rbind, lapply(names(exons.count), function(sample) {
  tryCatch(do.call('rbind', lapply(LOGS[grepl(sample, LOGS, fixed = T)], function(log) {
    file <- read.delim(file.path(FILE_PATH, 'logs', log), header = F, comment.char = '#', as.is = T)
    
    data.frame(total = file[5, 2] %>% as.integer,
               unique = file[8, 2] %>% as.integer,
               multimapped = file[25, 2] %>% as.integer,
               unmapped = c(file[28, 2], file[30, 2], file[32, 2]) %>% as.integer %>% sum) %>% as.matrix
  })) %>% colSums,
  error = function(message) {
    message(paste('No logs for', sample))
    data.frame(total = 0, unique = 0, multimapped = 0, unmapped = 0)
  })
})) %>% `rownames<-`(names(exons.count))

# QC ----------------------------------------------------------------------

qc <- qc %>% mutate(mitochondrial = mitochondrial / total, rRNA = rRNA / total) %>%
  `rownames<-`(rownames(qc)) %>% as.matrix %>% reshape2::melt(varnames = c('Sample', 'Metric')) %>%
  mutate(Metric = factor(Metric, levels = c('mitochondrial', 'rRNA', 'total'), labels = c('Mitochondrial', 'Ribosomal', 'Total')))

qc[qc$Metric != 'Total', ] %>%
  ggplot(aes(Sample, value, fill = Metric)) + 
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  ggtitle('Quality Metrics') + labs(tag = 'A') +
  theme(axis.text.x = element_blank()) +
  
  qc[qc$Metric == 'Total', ] %>%
  ggplot(aes(Sample, value, fill = Metric)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(strip.text = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('Total Reads') + labs(tag = 'B') +
  
  plot_layout(guides = 'collect', ncol = 1) &
  theme(legend.position = 'none', strip.background = element_blank()) &
  xlab(element_blank()) & ylab(element_blank()) & facet_wrap(~Metric, ncol = 1, scales = 'free_y') &
  geom_bar(stat = 'identity') -> grid

ggsave2(file.path(FILE_PATH, 'quantified', 'qc.pdf'), grid, width = 11, height = 6)

# Mapping -----------------------------------------------------------------

data.plot <- data.frame(Sample = names(exons.count),
                        Type = c(rep('Exonic', length(exons.count)), rep('Intronic', length(introns.count)),
                                 rep('Intergenic', length(exons.count)), rep('Multimapped', length(exons.count)),
                                 rep('Unmapped', length(exons.count))),
                        Reads = c(exons.count, introns.count,
                                  qc.STAR[, 'unique'] - exons.count - introns.count, qc.STAR[, 'multimapped'],
                                  qc.STAR[, 'unmapped'])) %>% mutate(Reads.fraction = Reads / qc.STAR[, 'total'])

data.plot %>%
  ggplot(aes(Sample, Reads.fraction, fill = Type)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  theme(axis.text.x = element_blank()) +
  ggtitle('Fraction') + labs(tag = 'A') +
  
  data.plot %>%
    ggplot(aes(Sample, Reads, fill = Type)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle('Raw') + labs(tag = 'B') +
  
  plot_layout(guides = 'collect', ncol = 1) &
  theme(legend.position = 'bottom') &
  xlab(element_blank()) & ylab(element_blank()) & geom_bar(stat = 'identity') -> grid

ggsave2(file.path(FILE_PATH, 'quantified', 'mapping_distribution.pdf'), grid, width = 11, height = 6)

# Introns/Exons -----------------------------------------------------------

data.plot <- data.frame(Sample = names(exons.count),
                        Location = c(rep('Exon', length(exons.count)), rep('Intron', length(introns.count))) %>% factor(c('Intron', 'Exon')),
                        Reads = c(exons.count, introns.count)) %>% mutate(Reads.fraction = Reads / (introns.count + exons.count))

data.plot %>%
  ggplot(aes(Sample, Reads.fraction, fill = Location)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  theme(axis.text.x = element_blank()) +
  ggtitle('Fraction') + labs(tag = 'A') +
  
  data.plot %>%
    ggplot(aes(Sample, Reads, fill = Location)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle('Raw') + labs(tag = 'B') +
  
  plot_layout(guides = 'collect', ncol = 1) &
  theme(legend.position = 'bottom') &
  xlab(element_blank()) & ylab(element_blank()) & geom_bar(stat = 'identity') -> grid

ggsave2(file.path(FILE_PATH, 'quantified', 'count_distribution.pdf'), grid, width = 11, height = 6)