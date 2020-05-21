# Initialize --------------------------------------------------------------

args <- commandArgs(T)

if(length(args) < 3) {
  stop('Usage: Rscript quantify.R <human|mouse> <paired> <input directory> [output directory] [cores]', call. = F)
} else if (!dir.exists(args[3])) {
  stop(paste0("'", args[3], "' is not a valid directory."), call. = F)
}

SPECIES <- args[1]

if(!(SPECIES %in% c('human', 'mouse')))
  stop('Usage: Rscript quantify.R <human|mouse> <paired> <input directory> [output directory] [cores]', call. = F)

PAIRED <- args[2] == 'true'

IN_DIR <- args[3]
OUT_DIR <- args[3]

if(length(args) > 3)
  OUT_DIR <- args[4]

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicAlignments)
  library(edgeR)
  library(rtracklayer)
  library(dplyr)
  library(Matrix)
  library(parallel)
})

CORES <- 5
if(length(args) > 4) {
  CORES <- suppressWarnings(as.integer(args[5]))
  if(is.na(CORES)) {
    message('Cores (argument 5) was non-integer. Defaulting to 5.')
  }
}
options(mc.cores = CORES)

if(!dir.exists(file.path(OUT_DIR, 'quantified', 'logs')))
  dir.create(file.path(OUT_DIR, 'quantified', 'logs'), recursive = T)

# Build the reference if it doesn't exist already.
if(!file.exists(file.path('../references', paste0(SPECIES, '.Rdata')))) {
  message('Reference didn\'t exist. Building.')
  system(paste('Rscript buildReference.R', SPECIES))
}

# Load the genome information
load(file.path('../references', paste0(SPECIES, '.Rdata')))
PARAM.BAM <- ScanBamParam(tag = 'NH', what = c('qname', 'flag'))

MITOCHONDRIAL <- EXONS[Filter(function(gene) gene %in% names(EXONS), ANNOTATIONS$gene_name[ANNOTATIONS$mitochondrial])]
RRNA <- EXONS[Filter(function(gene) gene %in% names(EXONS), ANNOTATIONS$gene_name[ANNOTATIONS$gene_biotype == 'rRNA'])]

# Process -----------------------------------------------------------------

FILES <- list.files(IN_DIR, '\\.bam$')

message(paste('Starting quantification of', length(FILES), 'files.'))

# Generate a list of sparse Matrices (for each sample) that has quantified introns,
# exons, etc. in columns and gene IDs on rows.
COUNTS <- mclapply(FILES, function(file) {
  # Return a tmp file if we have it.
  if(file.exists(file.path(OUT_DIR, 'quantified', paste0('tmp_', substring(file, 0, nchar(file) - 4), '.rds'))))
    return(readRDS(file.path(OUT_DIR, 'quantified', paste0('tmp_', substring(file, 0, nchar(file) - 4), '.rds'))))
  
  file.create(file.path(OUT_DIR, 'quantified', 'logs', paste0(substring(file, 0, nchar(file) - 4), '.out')))
  
  cat(paste(format(Sys.time(), '%b %d %H:%M:%S'), '.....', 'starting quantification'),
      file = file.path(OUT_DIR, 'quantified', 'logs', paste0(substring(file, 0, nchar(file) - 4), '.out')),
      sep = '\n')
  
  if(PAIRED) {
    reads <- tryCatch(readGAlignmentPairs(file.path(IN_DIR, file), param = PARAM.BAM),
                      error = function(message) {
                        cat(paste(format(Sys.time(), '%b %d %H:%M:%S'), '.....', 'no reads mapped in BAM'),
                            file = file.path(OUT_DIR, 'quantified', 'logs', paste0(substring(file, 0, nchar(file) - 4), '.out')),
                            sep = '\n',
                            append = T)
                        NULL
                      })
  } else {
    reads <- tryCatch(readGAlignments(file.path(IN_DIR, file), param = PARAM.BAM),
                      error = function(message) {
                        cat(paste(format(Sys.time(), '%b %d %H:%M:%S'), '.....', 'no reads mapped in BAM'),
                            file = file.path(OUT_DIR, 'quantified', 'logs', paste0(substring(file, 0, nchar(file) - 4), '.out')),
                            sep = '\n',
                            append = T)
                        NULL
                      })
  }
  
  
  if(is.null(reads)) {
    counts <- data.frame(exon = 0, intron = 0) %>% as.matrix
    fpkm <- data.frame(exon = 0, intron = 0) %>% as.matrix
    counts.mito <- 0
    counts.rRNA <- 0
    total <- 0
  } else {
    if(PAIRED) {
      total <- length(unique(values(reads %>% first)$qname))
      reads.unique <- reads[values(reads %>% first)$flag < 1024]
    } else {
      total <- length(unique(values(reads)$qname))
      reads.unique <- reads[values(reads)$flag < 1024]
    }
    
    rm(reads)
    
    junc <- junctions(reads.unique)
    no_junc <- elementNROWS(junc) == 0
    intron.overlap <- overlapsAny(reads.unique, INTRONS, minoverlap = 3, ignore.strand = T)
    
    counts <- assay(summarizeOverlaps(EXONS, reads.unique[!no_junc | !intron.overlap],
                                      'IntersectionNotEmpty',
                                      ignore.strand = T))
    counts <- cbind(counts, 0)
    colnames(counts) <- c('exon', 'intron')
    
    # Ensure similar ordering and dimensions
    introns <- assay(summarizeOverlaps(INTRONS, reads.unique[no_junc & intron.overlap],
                                       'IntersectionNotEmpty',
                                       ignore.strand = T))
    counts[rownames(introns), 'intron'] <- introns
    rm(introns)
    
    fpkm <- cbind(switch(all(counts[, 'exon'] == 0) + 1, rpkm(counts[, 'exon'], EXONS.SIZES), counts[, 'exon']),
                  switch(all(counts[, 'intron'] == 0) + 1, rpkm(counts[, 'intron'], INTRONS.SIZES), counts[, 'intron']))
    colnames(fpkm) <- c('exon', 'intron')
    
    if(PAIRED) {
      counts.mito <- length(unique(values(subsetByOverlaps(reads.unique, MITOCHONDRIAL, ignore.strand = T) %>% first)$qname))
      counts.rRNA <- length(unique(values(subsetByOverlaps(reads.unique, RRNA, ignore.strand = T) %>% first)$qname))
    } else {
      counts.mito <- length(unique(values(subsetByOverlaps(reads.unique, MITOCHONDRIAL, ignore.strand = T))$qname))
      counts.rRNA <- length(unique(values(subsetByOverlaps(reads.unique, RRNA, ignore.strand = T))$qname))
    }
  }
  
  cat(paste(format(Sys.time(), '%b %d %H:%M:%S'), '.....', 'done'),
      file = file.path(OUT_DIR, 'quantified', 'logs', paste0(substring(file, 0, nchar(file) - 4), '.out')),
      sep = '\n',
      append = T)
  
  ret <- list(counts = counts %>% Matrix(sparse = T),
              fpkm = fpkm %>% Matrix(sparse = T),
              mitochondrial = counts.mito, rRNA = counts.rRNA, total = total)
  saveRDS(ret, file.path(OUT_DIR, 'quantified', paste0('tmp_', substring(file, 0, nchar(file) - 4), '.rds')))
  ret
})

names(COUNTS) <- unlist(lapply(FILES, tools::file_path_sans_ext))

# Save --------------------------------------------------------------------

if(length(COUNTS) > 0) {
  # Save aggregated intron/exon as their corresponding file
  # Gene IDs will be on rows, samples will be on columns
  lapply(c('counts', 'fpkm'), function(type) {
    lapply(c('intron', 'exon'), function(pivot) {
      do.call(cbind, lapply(COUNTS, function(sample) {
        sample[[type]][, pivot]
      })) %>% as.matrix %>% Matrix(sparse = T) %>% `colnames<-`(names(COUNTS)[sample]) %>%
        saveRDS(file.path(OUT_DIR, 'quantified', paste0(type, '_', pivot, '.rds')))
    })
  })
  
  # Save quality metrics (number mapping to mitochondrial genes, rRNA and total reads)
  do.call(cbind, lapply(c('mitochondrial', 'rRNA', 'total'), function(pivot) {
    sapply(names(COUNTS), function(sample) tryCatch(COUNTS[[sample]][[pivot]], error = function(message) 0))
  })) %>% `colnames<-`(c('mitochondrial', 'rRNA', 'total')) %>% write.table(file.path(OUT_DIR, 'quantified', 'qc.tsv'), quote = F, sep = '\t')
}