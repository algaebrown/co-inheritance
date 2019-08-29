import pandas as pd

def best_e_value(target_genome_blast_result):
    '''
    :param target_genome_file:
    :return: best evalue of each gene within the target genome (series)
    '''
    df = pd.read_csv(target_genome_blast_result, names = ['query_gene', 'e_value'])

    # find best_evalue(smallest is the best)
    df.sort_values(by = 'e_value', axis = 'index', ascending = True, inplace = True)
    df.drop_duplicates(subset = 'query_gene', keep = 'first', inplace = True)

    # find smallest e_value among gene that is not zero (for log transform later)
    smallest_e_value = df.loc[df['e_value'] != 0]['e_value'].min() # the dataframe is sorted

    return(df, smallest_e_value)

from os import listdir
from os.path import isfile, join, exists
def append_to_file(all_results_folder, output_file):
    '''
    To concat all diamond blastp results in order to make pivot matrix
    To generate smallest non-zero e-value among all query-target pairs
    :param all_results_folder: folder containing all diamond blastp results
    :param output_file: concated blastp results with best evalue for each target genome
    :return: smallest_e_value
    '''
    # add header to output file
    with open(output_file, 'w') as o:
        o.write(','.join(['query_gene','e_value', 'target_genome'])+'\n')

    # initialize smallest_evalue with 1
    smallest_among_all = 1

    files = [f for f in listdir(all_results_folder) if isfile(join(all_results_folder, f))]
    for target_genome_blast_result in files:
        df, smallest_e_value = best_e_value(target_genome_blast_result)

        # add target genome information to df
        df['target_genome'] = target_genome_blast_result.split("\\")[-1] # retain target genome ID

        print(df.iloc[0])
        with open(output_file, 'a') as o:
            df.to_csv(o, header = False, index = False)

        # compete smallest e-value, update
        if smallest_among_all > smallest_e_value:
            smallest_among_all = smallest_e_value
    return (smallest_among_all)

import math
def log_e_value(e_value, smallest_among_all):
    '''
    Perform log transformation on e-value. See https://doi.org/10.1371/journal.pone.0139006
    :param e_value:
    :param smallest_among_all:
    :return: transformed e-value ranging from 0-1. 1 means perfect match. 0 means not match
    '''
    if e_value == 0: # perfect fit
        return(1)
    elif e_value >= 1: # screwed
        return(0)
    else:
        return (math.log(a) / -(smallest_among_all))

def transform_gene_wise(output_file, gene_wise_folder, bins = 100):
    '''
    To perform log transformation on all target-query pairs
    To organize into one gene per file
    To retrun global discretization
    :param output_file: The concatenated file with all blast results
    :param gene_wise_folder: The folder to contain one gene per file
    :return: global_qbins: quantile discretization
    '''

    # read concatenated file chunkwise:
    chunks = pd.read_csv(output_file, chunksize = 10000, header = 0)

    # initialize sampled trans_e_value
    sampled = []
    for chunk in chunks:

        # transfer e-value and discard original value
        chunk['trans_e_value'] = chunk['e_value'].apply(log_e_value)
        chunk.drop(labels = 'e_value', axis = columns, inplace = True)

        # sample 1% tran-evalue for estimation of discretization
        sampled = sampled + chunk['trans_e_value'].sample(frac = 0.01).tolist()

        # iterate over all rows, organize them gene wise
        for row in chunk.iterrows():
            gene = row['query_gene']
            # if not already a file
            if exists(gene_wise_folder + gene):
                with open(gene_wise_folder + gene, 'a') as f:
                    f.write(row['target_genome'] + ',{:4d}'.format(row['trans_e_value']) + '\n') # 4 digits is more than enough
            else:
                with open(gene_wise_folder + gene, 'w') as f:
                    f.write('target_genome, trans_e_value\n')
                    f.write(row['target_genome'] + ',{:4d}'.format(row['trans_e_value']) + '\n')

        # quantile discretization bins
        global_qbins = pd.qcut(sampled, bins-2, labels = range(1, bins -1), retbins = True)
        return(global_qbins)




def to_pivot(gene_wise_folder, pivot_file):
    '''

    :param gene_wise_folder:
    :param pivot_file:
    :return:
    '''
    files = [f for f in listdir(gene_wise_folder) if isfile(join(gene_wise_folder, f))]
    for f in files:
        gene = f.split('\\')[-1]
