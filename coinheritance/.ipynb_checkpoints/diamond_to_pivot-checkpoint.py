import pandas as pd

from index import gene_id, target_genome_id

def initialize(fasta_file, diamond_folder):
    global gid, tid
    gid = gene_id(fasta_file)  ### args
    tid = target_genome_id(diamond_folder)  ### args

def best_e_value(target_genome_blast_result):
    """
    select best e-value for each gene WITHIN each target genome
    :param target_genome_blast_result:
    :return: best evalue of each gene within the target genome (series)
    """
    df = pd.read_csv(target_genome_blast_result, names=['query_gene', 'e_value'], sep= '\t', header = None)

    # find best_evalue(smallest is the best)
    df.sort_values(by='e_value', axis='index', ascending=True, inplace=True)
    df.drop_duplicates(subset='query_gene', keep='first', inplace=True)

    # remap to shorter gene_id
    df['query_gene'] = df['query_gene'].map(gid)

    # find smallest e_value among gene that is not zero (for log transform later)
    smallest_e_value = df.loc[df['e_value'] != 0]['e_value'].min()  # the dataframe is sorted

    return df, smallest_e_value


from os import listdir
from os.path import isfile, join, exists


def concat_diamond_result(diamond_folder, concated_file):
    """
    To concat all diamond blastp results in order to make pivot matrix
    To generate smallest non-zero e-value among all query-target pairs
    :param diamond_folder: folder containing all diamond blastp results
    :param concated_file: concated blastp results with best evalue for each target genome
    :return: smallest_e_value
    """
    # add header to output file
    with open(concated_file, 'w') as o:
        o.write(','.join(['query_gene', 'e_value', 'target_genome']) + '\n')

    # initialize smallest_e_value with 1
    global smallest_among_all
    smallest_among_all = 1

    files = [f for f in listdir(diamond_folder) if isfile(join(diamond_folder, f))]
    for target_genome_blast_result in files:

        df, smallest_e_value = best_e_value(diamond_folder + target_genome_blast_result)

        # add target genome information to df
        df['target_genome'] = tid[target_genome_blast_result]  # retain target genome ID

        print(df.iloc[0])
        with open(concated_file, 'a') as o:
            df.to_csv(o, header=False, index=False)

        # compete smallest e-value, update
        if smallest_among_all > smallest_e_value:
            smallest_among_all = smallest_e_value

    return smallest_among_all


import math


def log_e_value(e_value):
    """
    Perform log transformation on e-value. See https://doi.org/10.1371/journal.pone.0139006
    :param e_value:
    :param smallest_among_all:
    :return: transformed e-value ranging from 0-1. 1 means perfect match. 0 means not match
    """
    if e_value == 0:  # perfect fit
        return 1
    elif e_value >= 1:  # screwed
        return 0
    else:
        return math.log(e_value) / (- math.log(smallest_among_all))

from scipy.sparse import csr_matrix, save_npz


import numpy as np
def to_sparse_matrix(concated_file):
    """

    :param smallest_among_all:
    :param concated_file:
    :return: scipy sparse matrix
    """


    row_index = []
    col_index = []
    data = []



    chunks = pd.read_csv(concated_file, chunksize=10000, header=0)

    for chunk in chunks:
        row_index = row_index + chunk['query_gene'].tolist()
        col_index = col_index + chunk['target_genome'].tolist()
        data = data + chunk['e_value'].apply(log_e_value).tolist()

    # initialize sparse matrix
    ma = csr_matrix((data, (row_index, col_index)), shape = (len(gid), len(tid)), dtype=np.half)  # row = gene, column = target_genome
    return ma, data

from random import sample
def qcut_bins(data, bins = 100):
    '''

    :param bins:
    :param data:
    :return:
    '''
    s = pd.qcut(sample(data, 10000), q = bins, labels = range(1, bins + 1), retbins = True)
    return s

def main(fasta_file, diamond_folder, tmp_folder, bins = 100):
    '''

    :param bins:
    :param fasta_file:
    :param diamond_folder:
    :param tmp_folder:
    :return:
    '''
    initialize(fasta_file, target_genome_id)
    concat_diamond_result(diamond_folder, tmp_folder + 'concated_diamond')
    ma, data = to_sparse_matrix(tmp_folder + 'concated_diamond')
    save_npz(tmp_folder + 'sparse_pivot.npz',ma, compressed = True)
    s = qcut_bins(data, bins = bins)
    np.save(tmp_folder + 'qcut_'+ str(bins)+ '.npz', s)
