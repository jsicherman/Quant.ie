library(Rsamtools)
library(GenomicAlignments)
library(edgeR)
library(rtracklayer)
library(dplyr)
library(Matrix)
library(parallel)

# Initialize --------------------------------------------------------------

args <- commandArgs(trailingOnly = T)

if(length(args) < 3) {
  stop('Usage: Rscript quantify.R <human|mouse> <paired> <input directory> [output directory] [cores]', call. = F)
} else if (!dir.exists(args[1])) {
  stop(paste0("'", args[1], "' is not a valid directory."), call. = F)
}

SPECIES <- args[1]

if(!(species %in% c('human', 'mouse')))
  stop('Usage: Rscript quantify.R <human|mouse> <paired> <input directory> [output directory] [cores]', call. = F)

IN_DIR <- args[2]
OUT_DIR <- args[2]

PAIRED <- args[3] == 'true'

if(length(args) > 3)
  OUT_DIR <- args[4]

if(length(args) > 4) {
  CORES <- suppressWarnings(as.integer(args[5]))
  if(is.na(CORES)) {
    message('Cores (argument 5) was non-integer. Defaulting to 5.')
    CORES <- 5
  }
}
options(mc.cores = CORES)

if(!dir.exists(OUT_DIR))
  dir.create(OUT_DIR, recursive = T)

# Build the reference if it doesn't exist already.
if(!file.exists(file.path('references', paste0(SPECIES, '.Rdata'))))
  system(paste('Rscript src/buildReference.R', SPECIES))

# Load the genome information
load(file.path('references', paste0(SPECIES, '.Rdata')))
PARAM.BAM <- ScanBamParam(tag = 'NH', what = c('qname', 'flag'))

MITOCHONDRIAL <- EXONS[Filter(function(gene) gene %in% names(EXONS), ANNOTATIONS$gene_id[ANNOTATIONS$mitochondrial])]
RRNA <- EXONS[Filter(function(gene) gene %in% names(EXONS), ANNOTATIONS$gene_id[ANNOTATIONS$gene_biotype == 'rRNA'])]

# Process -----------------------------------------------------------------

FILES <- list.files(IN_DIR, '\\.bam$')

# Generate a list of sparse Matrices (for each sample) that has quantified introns,
# exons, etc. in columns and gene IDs on rows.
COUNTS <- mclapply(FILES, function(file) {
  reads <- tryCatch(ifelse(PAIRED, readGAlignmentsPairs, readGAlignments)(file.path(IN_DIR, file), param = PARAM.BAM),
                    error = function(message) {
                      message('There are no genome reads mapped in BAM') 
                      NULL
                    })
  
  if(is.null(reads)) {
    counts.exon <- 0
    counts.intron <- 0
    counts.mito <- 0
    counts.rRNA <- 0
  } else {
    if(PAIRED) {
      total <- length(unique(values(reads %>% first)$qname))
      reads.unique <- reads[values(reads %>% first)$flag < 1024]
    } else {
      total <- length(unique(values(reads)$qname))
      reads.unique <- reads[values(reads)$flag < 1024]
    }
    
    junc <- junctions(reads.unique)
    no_junc <- elementNROWS(junc) == 0
    intron.overlap <- overlapsAny(reads.unique, INTRONS, minoverlap = 3, ignore.strand = T)
    
    counts.intron <- assay(summarizeOverlaps(INTRONS, reads.unique[no_junc & intron.overlap],
                                                'IntersectionNotEmpty',
                                                ignore.strand = T))
    
    counts.exon <- assay(summarizeOverlaps(EXONS, reads.unique[!no_junc | !intron.overlap],
                                          'IntersectionNotEmpty',
                                          ignore.strand = T))
    
    counts.mito <- length(unique(values(subsetByOverlaps(reads, MITOCHONDRIAL, ignore.strand = T))$qname))
    counts.rRNA <- length(unique(values(subsetByOverlaps(reads, RRNA, ignore.strand = T))$qname))
  }
  
  list(introns = counts.intron, exons = counts.exon,
       mitochondrial = counts.mito, rRNA = counts.rRNA,
       total = total)
})

names(COUNTS) <- unlist(lapply(FILES, tools::file_path_sans_ext))

# Save --------------------------------------------------------------------

if(length(COUNTS) > 0) {
  # Save aggregated intron/exon as their corresponding file
  # Gene IDs will be on rows, samples will be on columns
  lapply(c('introns', 'exons'), function(pivot) {
    do.call(cbind, lapply(names(COUNTS), function(sample)
      COUNTS[[sample]][[pivot]]
    )) %>% as.matrix %>% Matrix(sparse = T) %>%
      `colnames<-`(names(COUNTS)) %>%
      saveRDS(file.path(OUT_DIR, paste0(pivot, '.rds')))
  })
  
  # Save quality metrics (number mapping to mitochondrial genes, rRNA and total reads)
  do.call(cbind, lapply(c('mitochondrial', 'rRNA', 'total'), function(pivot) {
    sapply(names(COUNTS), function(sample) COUNTS[[sample]][[pivot]])
  })) %>% `colnames<-`(c('mitochondrial', 'rRNA', 'total')) %>% write.table(file.path(OUT_DIR, 'qc.tsv'), quote = F, sep = '\t')
}