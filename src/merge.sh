#!/bin/bash

IN_DIR=""
OUT_DIR=""
N_THREADS=10

while getopts ":i:o:ht:" opt; do
  case $opt in
    i) IN_DIR="$OPTARG"
    ;;
    o) OUT_DIR="$OPTARG"
    ;;
    t) N_THREAD="$OPTARG"
    ;;
    h) echo "usage: merge -i input-directory -o output-directory"
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

mkdir -p $OUT_DIR

OUT_DIR="$(realpath ${OUT_DIR})"
IN_DIR="$(realpath ${IN_DIR})"
cd $IN_DIR

# Get all samples, identified by their first lane
files=()
for file in *_L001.bam; do
  [ -f "$file" ] || break
  
  files+=(`echo "$file" | rev | cut -d_ -f2- | rev`)
done

prog() {
  local w=80 p=$1 t=$2 shift
  printf -v dots "%*s" "$(( $p*$w/$t ))" ""; dots=${dots// /#};
  printf "\r\e[K|%-*s| %3d %% %s" "$w" "$dots" "$(($p*100/$t))"; 
}

cd $OUT_DIR

# Merge BAM files for same samples on different lanes
echo "Merging BAM files for ${#files[@]} samples."
for ((i=0;i<${#files[@]};++i)); do
  prog "$i" "${#files[@]}"
  
  samtools merge -@ $N_THREADS "${files[$i]}.bam" `find "$IN_DIR" -maxdepth 1 -type f -name "${files[$i]}_*"`
done

prog 1 1
printf "\n"