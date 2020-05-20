IN_DIR=""
OUT_DIR=""

while getopts ":g:i:o:pht:n" opt; do
case $opt in
  i) IN_DIR="$OPTARG"
  ;;
  o) OUT_DIR="$OPTARG"
  ;;
  h) echo "usage: qc -i input-directory -o output-directory"
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

mkdir -p "$OUT_DIR"

OUT_DIR="$(realpath $OUT_DIR)"
cd $IN_DIR

rm -f $OUT_DIR/fastqc.sh

echo "#!/bin/bash" >> $OUT_DIR/fastqc.sh
echo cd "$(realpath ./)" >> $OUT_DIR/fastqc.sh

batch=1
i=0
for file in *.fastq*; do
  [ -f "$file" ] || break
  
  if [ "$(( $i % 10 ))" -eq "0" ]; then
    echo "if [ \$SLURM_ARRAY_TASK_ID = $batch ]; then" >> $OUT_DIR/fastqc.sh
  fi
  
  echo "  fastqc -o $OUT_DIR $file" >> $OUT_DIR/fastqc.sh
  
  if [ "$(( ($i + 1) % 10 ))" -eq "0" ]; then
    echo "fi" >> $OUT_DIR/fastqc.sh
    ((batch++))
  fi
  
  ((i++))
done

if [ "$(( $i % 10 ))" -ne "0" ]; then
  echo "fi" >> $OUT_DIR/fastqc.sh
  ((batch++))
fi

echo "if [ \$SLURM_ARRAY_TASK_ID = $batch ]; then" >> $OUT_DIR/fastqc.sh
echo "  multiqc ." >> $OUT_DIR/fastqc.sh
echo "fi" >> $OUT_DIR/fastqc.sh

sbatch --array "1-$batch" $OUT_DIR/fastqc.sh