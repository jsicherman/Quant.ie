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

echo "#!/bin/bash" >> $OUT_DIR/fastqc.sh
echo cd "$(realpath ./)" >> $OUT_DIR/fastqc.sh

i=1
for file in *.fastq*; do
  [ -f "$file" ] || break
  
  echo "if [ \$SLURM_ARRAY_TASK_ID = $i ]; then" >> $OUT_DIR/fastqc.sh
  echo "  fastqc -o $OUT_DIR $file" >> $OUT_DIR/fastqc.sh
  echo "fi" >> $OUT_DIR/fastqc.sh
  ((i++))
done

echo "if [ \$SLURM_ARRAY_TASK_ID = $i ]; then" >> $OUT_DIR/fastqc.sh
echo "  multiqc ." >> $OUT_DIR/fastqc.sh
echo "fi" >> $OUT_DIR/fastqc.sh

sbatch --array "1-$i" $OUT_DIR/fastqc.sh