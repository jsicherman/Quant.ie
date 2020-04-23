#!/bin/bash

GENOME="human"
GENOME_DIR="/cosmos/data/pipeline-output/rnaseq/references/mm10_ncbi"
IN_DIR=""
OUT_DIR=""
PAIRED=false
N_THREAD=5

while getopts ":g:i:o:pht:" opt; do
  case $opt in
    g) GENOME="$OPTARG"
    case $GENOME in
      human) GENOME_DIR="/cosmos/data/pipeline-output/rnaseq/references/mm10_ncbi"
      ;;
      mouse) GENOME_DIR="/cosmos/data/pipeline-output/rnaseq/references/hg38_ncbi"
      ;;
    esac
    ;;
    i) IN_DIR="$OPTARG"
    ;;
    o) OUT_DIR="$OPTARG"
    ;;
    p) PAIRED=true
    ;;
    t) N_THREAD="$OPTARG"
    ;;
    h) echo "usage: quantify [-g [human|mouse]] [-p] [-t threads] -i input-directory -o output-directory"
    exit
    ;;
    \?) echo "Invalid option -$OPTARG" 1>&2
    exit 2
    ;;
  esac
done

if [ -z "$IN_DIR" ] || [ -z "$OUT_DIR" ]; then
  echo "Input/output directories must be specified" 1>&2
  exit 3
fi

if [ "$GENOME" != "human" ] && [ "$GENOME" != "mouse" ]; then
  echo "-g must be one of human or mouse" 1>&2
  exit 4
fi

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"
rm tasks.sh

if [ "$PAIRED" = true ]; then
  echo "Running paired end alignment"
else
  for file in $(find $IN_DIR -maxdepth 1 -type f | sort); do
    fileQualified="${file##*/}"
    fileName="${fileQualified%.*}"
    
    rm "$fileName.sh"
    
    echo "#!/bin/bash" >> $fileName.sh
    echo source /cvmfs/soft.computecanada.ca/config/profile/bash.sh >> $fileName.sh
    
    echo cd "$(realpath $OUT_DIR)" >> $fileName.sh
    echo module load star/2.7.3a >> $fileName.sh
    
    echo STAR --genomeDir "$GENOME_DIR" --genomeLoad LoadAndKeep --runThreadN $N_THREAD --readFilesIn $file --outFileNamePrefix "$fileName/" --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --quantMode GeneCounts --outFilterMultimapNmax 1 >> $fileName.sh
    echo STAR --inputBAMfile "$fileName/Aligned.sortedByCoord.out.bam" --runThreadN $N_THREAD --bamRemoveDuplicatesType UniqueIdentical --runMode inputAlignmentsFromBAM --outFileNamePrefix "$fileName/" >> $fileName.sh
    
    chmod +x $fileName.sh
  done
fi

echo "#!/bin/bash" >> tasks.sh
echo cd "$(realpath $OUT_DIR)" >> tasks.sh
i=1
for file in *.sh; do
  [ -f "$file" ] || break
  
  echo "if [ \$SLURM_ARRAY_TASK_ID = $i ]; then" >> tasks.sh
  echo "  ./$file" >> tasks.sh
  echo "fi" >> tasks.sh
  ((i++))
done

echo "Processed $(($i-1)) files. Running STAR."

sbatch --array "1-$i%2" tasks.sh