# Quant.ie
Pipeline for the quantification of introns and exons. This repository is a reproduction of [PSQ Pipeline](https://github.com/sonnyc247/PSQ_Pipeline/) written by Sonny Chen to work on the Pavlidis lab servers, as well as incorporating some optimizations and tweaks.

## Usage
The pipeline can be run in four steps using the following commands. All assume you're executing commands from within the src directory on your machine.

### Alignment using STAR
`./align.sh [-g [human|mouse]] [-p] [-t threads] -i input-directory -o output-directory`

| Flag | Description | Default |
| ---- | ----------- | ------- |
| g | Select the genome to align to. The specific reference genome used is either hg38_ensembl (human) or mm10_ensembl (mouse). | human |
| p | Whether or not to samples are paired-end. | false |
| t | How many threads to run STAR with. | 10 |
| i | The path where your input (fastq.gz) files are. | |
| o | The path where your output files will be sent. Will be created if it doesn't exist. | |

Files will be output to the following directories:

/output-directory/
|-- tasks.sh
|-- logs
|-- scripts
`-- bam
    |-- raw
    `-- processed

Each folder will contain one file per fastq.gz input (either a STAR log file, quantification shell script, bam file (without removing PCR duplicates) and bam file (after removing PCR duplicates), respectively).

If you run into problems trying to execute `./align.sh`, try running `chmod +x align.sh` to ensure you have execute permission on the file.

### Merging BAM
`./merge.sh [-t threads] -i input-directory -o output-directory`

| Flag | Description | Default |
| ---- | ----------- | ------- |
| t | How many additional threads to use. | 10 |
| i | The path where your input (.bam) files are. | |
| o | The path where your output files will be sent. Will be created if it doesn't exist. | |

BAM files for the same sample (along multiple lanes) will be merged using samtools. If you run into problems trying to execute `./merge.sh`, try running `chmod +x merge.sh` to ensure you have execute permission on the file.

### Quantification
`Rscript quantify.R <human|mouse> <paired> <input directory> [output directory] [cores]`

| Parameter | Description | Default |
| --------- | ----------- | ------- |
| 1 | What taxon you're interested in | |
| 2 | Whether or not your sample was paired | |
| 3 | The directory containing your BAM files | |
| 4 | The path to where you want your output. Should usually be the same as the output directory from `align.sh`. | Same as input directory |
| 5 | How many cores to use while doing quantification | 5 |

### Visualization
`Rscript visualize.R <directory>`

Echoes the percentage of reads that aligned to mitochondrial and ribosomal genes and saves figures of count distributions to quantified directory. The directory provided should be the parent folder to the count matrices. Your directory hierarchy should ideally look like the following.

/output-directory/
|-- tasks.sh
|-- logs
|-- scripts
|-- bam
    |-- raw
    `-- processed
|-- bam_merged
    |-- raw
    `-- processed
`quantified
    |-- counts_exon.rds
    |-- counts_intron.rds
    |-- fpkm_exon.rds
    `-- fpkm_intron.rds