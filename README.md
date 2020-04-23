# Quant.ie
Pipeline for the quantification of introns and exons.

## Usage
./quantify.sh [-g [human|mouse]] [-p] [-t threads] -i input-directory -o output-directory

| Flag | Description | Default |
| ---- | ----------- | ------- |
| g | Select the genome to align to. The specific reference genome used is either hg38_ncbi (human) or mm10_ncbi (mouse). | human |
| p | Whether or not to samples are paired-end. | false |
| t | How many threads to run STAR with. | 5 |
| i | The path where your input (fastq.gz) files are. | |
| o | The path where your output files will be sent. Will be created if it doesn't exist. | |