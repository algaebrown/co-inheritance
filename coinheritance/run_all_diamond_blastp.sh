#!/bin/bash
# execute in the folder containing all target genome .faa
# specify argument $1 = folder containing all target genome $2 = output folder $3 cdhit repr gene

# example usage: bash run_all_diamond_blastp.sh TARGET_GENOME_FOLDER OUTPUT_FOLDER REPRESENTING_GENES.faa

cd $1
pwd

for file in *.faa
do diamond makedb --in $file --db ${file%.faa}
done

for file in *.dmnd
do diamond blastp --query $3 --db ${file%.dmnd} --out $2${file%.dmnd} --outfmt 6 qseqid evalue --verbose --threads 8
done
