# Quant.ie
Pipeline for the quantification of introns and exons. This repository is a reproduction of [PSQ Pipeline](https://github.com/sonnyc247/PSQ_Pipeline/) written by Sonny Chen to work on the Pavlidis lab servers, as well as incorporating some optimizations and tweaks.

## Usage

### Alignment using STAR
`./align.sh [-g [human|mouse]] [-p] [-t threads] -i input-directory -o output-directory`

| Flag | Description | Default |
| ---- | ----------- | ------- |
| g | Select the genome to align to. The specific reference genome used is either hg38_ncbi (human) or mm10_ncbi (mouse). | human |
| p | Whether or not to samples are paired-end. | false |
| t | How many threads to run STAR with. | 5 |
| i | The path where your input (fastq.gz) files are. | |
| o | The path where your output files will be sent. Will be created if it doesn't exist. | |

### Quantification
`Rscript quantify.R <human|mouse> <paired> <input directory> [output directory] [cores]`

| Parameter | Description | Default |
| --------- | ----------- | ------- |
| 1 | What taxon you're interested in | |
| 2 | Whether or not your sample was paired | |
| 3 | The directory containing your BAM files | |
| 4 | The path to where you want your output | Same as input directory |
| 5 | How many cores to use while doing quantification | 5 |