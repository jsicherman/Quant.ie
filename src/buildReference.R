suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(rtracklayer)
  library(dplyr)
  library(data.table)
})

# Rscript buildReference.R <human|mouse>

if(length(commandArgs(T)) == 0)
  stop('Usage: Rscript buildReference.R <human|mouse>')

SPECIES <- commandArgs(T)[1]

if(!(SPECIES %in% c('human', 'mouse')))
  stop('Usage: Rscript buildReference.R <human|mouse>')

REFERENCE <- import.gff(ifelse(SPECIES == 'mouse',
                               '/space/grp/Pipelines/rnaseq-pipeline/Assemblies/mm10_ensembl98/Mus_musculus.GRCm38.98.gtf',
                               '/space/grp/Pipelines/rnaseq-pipeline/Assemblies/hg38_ensembl98/Homo_sapiens.GRCh38.98.gtf'))

# Keep all chromosomes (1-20, X, MT), ie. dropping contigs/patches
# NCBI: REFERENCE <- keepSeqlevels(REFERENCE, grep('^NC_.*', levels(seqnames(REFERENCE)), value = T), pruning.mode = 'coarse')
REFERENCE <- keepSeqlevels(REFERENCE, levels(seqnames(REFERENCE))[1:22], pruning.mode = 'coarse')

# Subset our data for exons/genes
REFERENCE.EXONS <- REFERENCE[REFERENCE$type == 'exon']
REFERENCE.GENES <- REFERENCE[REFERENCE$type == 'gene']

rm(REFERENCE)

# Extend 3' UTR -----------------------------------------------------------

FLANK <- 100
exons <- split(REFERENCE.EXONS, REFERENCE.EXONS$transcript_id)

# An efficient method for finding the extreme endings on the positive strand
P3UTR <- strand(REFERENCE.EXONS) == '+' & end(REFERENCE.EXONS) ==
  merge(data.table(transcript_id = REFERENCE.EXONS$transcript_id),
        data.table(transcript_id = names(exons),
                   end = sapply(end(exons), max, na.rm = T)),
        by = 'transcript_id', sort = F)$end

# Copied for the negative sense strand
N3UTR <- strand(REFERENCE.EXONS) == '-' & start(REFERENCE.EXONS) ==
  merge(data.table(transcript_id = REFERENCE.EXONS$transcript_id),
        data.table(transcript_id = names(exons),
                   start = sapply(start(exons), min, na.rm = T)),
        by = 'transcript_id', sort = F)$start

# Extend
end(REFERENCE.EXONS[P3UTR]) <- end(REFERENCE.EXONS[P3UTR]) + FLANK
start(REFERENCE.EXONS[N3UTR]) <- start(REFERENCE.EXONS[N3UTR]) - FLANK

rm(P3UTR, N3UTR, exons, FLANK)

# Gene Info ---------------------------------------------------------------

# Annotations
# NCBI: ANNOTATIONS <- REFERENCE.GENES[, c('gene_id', 'description', 'gene_biotype', 'pseudo')] %>% as.data.frame
# NCBI: ANNOTATIONS$pseudo <- ANNOTATIONS$pseudo %>% factor(c(NA, 'true'), c('F', 'T'), exclude = NULL)
ANNOTATIONS <- REFERENCE.GENES[, c('gene_id', 'gene_name', 'gene_biotype')] %>% as.data.frame
ANNOTATIONS$mitochondrial <- ANNOTATIONS$seqnames == 'MT'
ANNOTATIONS$gene_biotype <- ANNOTATIONS$gene_biotype %>% as.factor
ANNOTATIONS[, c('seqnames', 'start', 'end', 'width')] <- NULL

# Exons
EXONS <- reduce(split(REFERENCE.EXONS, REFERENCE.EXONS$gene_name)) %>% as('GRangesList')
EXONS.SIZES <- sum(width(EXONS))

# Introns, mapped to overlapping genes
strand(REFERENCE.GENES) <- '*'
strand(REFERENCE.EXONS) <- '*'
INTRONS <- setdiff(REFERENCE.GENES, REFERENCE.EXONS)
intron.map <- findOverlaps(INTRONS, REFERENCE.GENES, type = 'within', select = 'all') %>% as.matrix
INTRONS <- reduce(split(INTRONS[intron.map[, 1]], REFERENCE.GENES$gene_name[intron.map[, 2]])) %>% as('GRangesList')
rm(intron.map)

INTRONS.SIZES <- EXONS.SIZES
INTRONS.SIZES[names(INTRONS)] <- sum(width(INTRONS))

# And save
if(!dir.exists('../references'))
  dir.create('../references', recursive = T)
save(INTRONS, EXONS, INTRONS.SIZES, EXONS.SIZES, ANNOTATIONS, file = file.path('../references', paste0(SPECIES, '.Rdata')))