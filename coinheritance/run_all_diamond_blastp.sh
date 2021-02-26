#!/bin/bash
# execute in the folder containing all target genome .faa
# specify argument $1 = folder containing all target genome $2 cdhit repr gene $3 outdir

# example usage: bash run_all_diamond_blastp.sh TARGET_GENOME_FOLDER OUTPUT_FOLDER REPRESENTING_GENES.faa
target_genome_dir=$1
gene_faa=$2
outdir=$3
run_diamond(){
    # Compare representing gene to all reference genome using diamond
    cd $target_genome_dir
    pwd
    # build diamond database
    for file in *.faa.gz
    do diamond makedb --in $file --db ${file%.faa.gz}
    done
    # make a new dir
    mkdir $outdir/diamond_results

    # compare
    for file in *.dmnd
    do diamond blastp --query $gene_faa --db ${file%.dmnd} --out $outdir/diamond_results/${file%.dmnd} --outfmt 6 qseqid evalue --verbose --threads 8
    done
    #remove all databases
    rm *.dmnd

}
run_diamond()
#python diamond_to_pivot.py -r $gene_faa -d $outdir/diamond_results -o $outdir
#python cal_score.py --table $outdir/binned.pivot.npz -o $outdir