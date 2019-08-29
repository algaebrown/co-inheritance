#!/bin/bashx
# $1=output file
# $2=folder containing blast result

# example usage: bash concat_blastp.sh CONCATED_FILE FOLDER_WITH_BLAST_RESULT

cd $2
for file in *
do
    awk '{print FILENAME (NF?"\t":"") $0}' $file >> $1
done
