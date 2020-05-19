#!/bin/bash

STAR_PATH="/space/grp/Pipelines/rnaseq-pipeline/Requirements/STAR/bin/Linux_x86_64/"
GENOME="human"
GENOME_DIR="/cosmos/data/pipeline-output/rnaseq/references/mm10_ensembl98"
IN_DIR=""
OUT_DIR=""
PAIRED=false
N_THREAD=10

while getopts ":g:i:o:pht:" opt; do
  case $opt in
    g) GENOME="$OPTARG"
    case $GENOME in
      human) GENOME_DIR="/cosmos/data/pipeline-output/rnaseq/references/mm10_ensembl98"
      ;;
      mouse) GENOME_DIR="/cosmos/data/pipeline-output/rnaseq/references/hg38_ensembl98"
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
    h) echo "usage: align [-g [human|mouse]] [-p] [-t threads] -i input-directory -o output-directory"
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

mkdir -p "$OUT_DIR"/scripts
mkdir -p "$OUT_DIR"/logs
mkdir -p "$OUT_DIR"/bam/raw
mkdir -p "$OUT_DIR"/bam/processed
cd "$OUT_DIR"
rm -f tasks.sh

prog() {
  local w=80 p=$1 t=$2 shift
  printf -v dots "%*s" "$(( $p*$w/$t ))" ""; dots=${dots// /#};
  printf "\r\e[K|%-*s| %3d %% %s" "$w" "$dots" "$(($p*100/$t))"; 
}

if [ "$PAIRED" = true ]; then
  echo "Generating scripts for paired end alignment"

  # Get each direction in a different array
  forward=()
  for file in `find "$IN_DIR" -maxdepth 1 -type f | sort | awk 'NR % 2 == 1'`; do
    forward+=("$file")
  done

  reverse=()
  for file in `find "$IN_DIR" -maxdepth 1 -type f | sort | awk 'NR % 2 == 0'`; do
    reverse+=("$file")
  done

  for ((i=0;i<${#forward[@]};++i)); do
    prog "$i" "${#forward[@]}"

    fileQualifiedForward="${forward[$i]##*/}"
    fileNameForward="${fileQualifiedForward%.*}"
    
    fileQualifiedReverse="${reverse[$i]##*/}"
    fileNameReverse="${fileQualifiedReverse%.*}"
    
    fileName=`echo "$fileNameForward" | rev | cut -d_ -f3- | rev`
    
    # LoadAndRemove for the last three alignments. Should clean up shared memory, but should double check
    GENOME_LOAD=LoadAndKeep
    if [ "$i" -ge "$(( ${#forward[@]} - 3 ))" ]; then
      GENOME_LOAD=LoadAndRemove
    fi
    
    rm -f "scripts/$fileName.sh"
    
    echo "#!/bin/bash" >> scripts/$fileName.sh
    echo mkdir "$(realpath ./)/$fileName" >> scripts/$fileName.sh
    echo cd "$(realpath ./)/$fileName" >> scripts/$fileName.sh
    
    echo "$STAR_PATH"STAR --genomeDir $GENOME_DIR --genomeLoad $GENOME_LOAD --runThreadN $N_THREAD --readFilesIn "${forward[$i]}" "${reverse[$i]}" --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --limitBAMsortRAM 10000000000 --outFilterMultimapNmax 1 >> scripts/$fileName.sh
    echo "$STAR_PATH"STAR --inputBAMfile Aligned.sortedByCoord.out.bam --runThreadN $N_THREAD --bamRemoveDuplicatesType UniqueIdentical --runMode inputAlignmentsFromBAM >> scripts/$fileName.sh
    
    echo mv Log.final.out "../logs/$fileName.out" >> scripts/$fileName.sh
    echo mv Aligned.sortedByCoord.out.bam "../bam/raw/$fileName.bam" >> scripts/$fileName.sh
    echo mv Processed.out.bam "../bam/processed/$fileName.bam" >> scripts/$fileName.sh
    echo cd ../ >> scripts/$fileName.sh
    echo rm -r $fileName >> scripts/$fileName.sh
    
    chmod +x scripts/$fileName.sh
  done
  
  prog 1 1
  printf "\n"
else
  echo "Generating scripts for unpaired alignment"
  
  files=$(find $IN_DIR -maxdepth 1 -type f | sort)
  for file in $files; do
    prog "$i" "${#files[@]}"
    
    fileQualified="${file##*/}"
    fileName="${fileQualified%.*}"
    
    # Skip mate-pairs by file name in unpaired mode
    if [[ $fileName == *"_L002_"* ]]; then
      echo "Skipping mate-pair '${fileName}'"
      continue
    fi
    
    # LoadAndRemove for the last three alignments. Should clean up shared memory, but should double check
    GENOME_LOAD=LoadAndKeep
    if [ "$i" -ge "$(( ${#files[@]} - 2 ))" ]; then
      GENOME_LOAD=LoadAndRemove
    fi
    
    rm -f "scripts/$fileName.sh"
    
    echo "#!/bin/bash" >> scripts/$fileName.sh
    echo mkdir "$(realpath ./)/$fileName" >> scripts/$fileName.sh
    echo cd "$(realpath ./)/$fileName" >> scripts/$fileName.sh
    
    echo "$STAR_PATH"STAR --genomeDir $GENOME_DIR --genomeLoad $GENOME_LOAD --runThreadN $N_THREAD --readFilesIn $fileName --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --limitBAMsortRAM 10000000000 --outFilterMultimapNmax 1 >> scripts/$fileName.sh
    echo "$STAR_PATH"STAR --inputBAMfile Aligned.sortedByCoord.out.bam --runThreadN $N_THREAD --bamRemoveDuplicatesType UniqueIdentical --runMode inputAlignmentsFromBAM >> scripts/$fileName.sh
    
    echo mv Log.final.out "../logs/$fileName.out" >> scripts/$fileName.sh
    echo mv Aligned.sortedByCoord.out.bam "../bam/raw/$fileName.bam" >> scripts/$fileName.sh
    echo mv Processed.out.bam "../bam/processed/$fileName.bam" >> scripts/$fileName.sh
    echo cd ../ >> scripts/$fileName.sh
    echo rm -r $fileName >> scripts/$fileName.sh
    
    chmod +x $fileName.sh
  done
  
  prog 1 1
  printf "\n"
fi

echo "#!/bin/bash" >> tasks.sh
echo cd "$(realpath ./)" >> tasks.sh

echo "Building task master..."

i=1
for file in scripts/*.sh; do
  [ -f "$file" ] || break
  
  echo "if [ \$SLURM_ARRAY_TASK_ID = $i ]; then" >> tasks.sh
  echo "  ./$file" >> tasks.sh
  echo "fi" >> tasks.sh
  ((i++))
done

echo "Processed $(($i-1)) files. Running STAR."

# Allocate three Slurm nodes and run jobs spread across them
salloc --spread-job -w cn[06-08] sbatch --array "1-$i%3" tasks.sh