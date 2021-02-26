import pandas as pd
import os
import math
import numpy as np
from scipy.sparse import csr_matrix, lil_matrix, save_npz
from index import gene_id, target_genome_id

def read_diamond(target_genome_blast_result):
    ''' read diamond output file to pandas dataframe
    '''
    df = pd.read_csv(target_genome_blast_result, names = ['qseqid','evalue'], delimiter= '\t', header = None, index_col = None)

    # find best_evalue(smallest is the best)
    df.sort_values(by='evalue', axis='index', ascending=True, inplace=True)
    df.drop_duplicates(subset='qseqid', keep='first', inplace=True)
    return df

def best_e_value(df):
    """
    select best e-value for each gene WITHIN each target genome
    :param target_genome_blast_result:
    :return: best evalue of each gene within the target genome (series)
    """
    
    # find smallest e_value among gene that is not zero (for log transform later)
    smallest_e_value = df.loc[df['evalue'] != 0]['evalue'].min()  # the dataframe is sorted

    return smallest_e_value
def diamond_to_pivot(target_path, tid, gid):
    '''
    extract evalue, log normalize it, and save in scipy sparse matrix
    target_path: str, path containing all diamond blastp outputs
    tid: dict, mapping target name to integers (column in pivot table)
    gid: dict, mapping gene name to integer (row in pivot table)
    returns: np.lil_matrix (genes * genome)
    '''
    
    # initialize sparse matrix
    ma = lil_matrix((len(gid), len(tid)))
    # initialize normalization factor (second best evalue)
    sec_best_eval = 1
    
    target_files = os.listdir(target_path)
    for file in target_files:
        df = read_diamond(os.path.join(target_path, file))
        # change to gene id
        df['gid'] = df['qseqid'].map(gid)
        df.dropna(subset = ['gid'], inplace = True)
        df['gid'] = df['gid'].astype(int)
        
        # target id
        target_id = tid[file]
        
        log_eval = -np.log(df['evalue']) #np.inf means evalue is 0
        ma[df['gid'].tolist(), target_id] = log_eval.tolist() # store negative log to prevent underflow :)
        
        # best_evl
        local_best_eval = best_e_value(df)
        if local_best_eval < sec_best_eval:
            sec_best_eval = local_best_eval
        
    # normalize by sec best
    ma = ma/(-np.log(sec_best_eval))
        
    # replace inf with 1
    npm = ma.toarray()
    npm[np.where(npm == np.inf)] = 1
        
    ma = lil_matrix(npm)
        
    return ma, sec_best_eval
def discretize(ma, bins = 100):
    # convert to csr matrix and qcut
    csr = ma.tocsr()
    q = pd.qcut(csr.data, q = bins, labels = False, duplicates = 'drop')
    
    # create a new csr matrix
    binned = csr_matrix((q, csr.indices, csr.indptr))
    
    return binned

    


def main(fasta_file, diamond_folder, bins = 100, outdir = '/tmp'):
    '''

    :param bins:
    :param fasta_file:
    :param diamond_folder:
    :param tmp_folder:
    :return:
    '''
    # map genome/genes to id 
    print('Getting index file to genome and genes')
    gid = gene_id(fasta_file, outdir = outdir)
    tid = target_genome_id(diamond_folder, outdir = outdir)

    # get normalized e-value
    print('computing pivot table from DIAMOND outputs')
    m, e = diamond_to_pivot(diamond_folder, tid, gid)
    print('Output pivot shape {}'.format(str(m.shape)))

    # discretize into 100 bins
    print('binning into {} bins'.format(bins))
    binned = discretize(m, bins = bins)

    return m, binned

def option_parser():
    from optparse import OptionParser
    usage = """
        THIS IS CO-INHERITANCE 1.0.0
        python diamond_to_pivot.py -r <query_genes.faa> -d <folder> -o <outdir>
        python diamond_to_pivot.py -r /home/hermuba/data0118/cdhit/Escherichia0.70rm_plasmid -d /home/hermuba/data0118/phylo_profile_archaea -o ~/data0118 """
    description = """CO-INHERTITANCE is a tool to compute co-inheritance relationship for a set of genes.
                    methods is by this paper:
                    Shin J, Lee I. Co-Inheritance Analysis within the Domains of Life Substantially Improves Network Inference by Phylogenetic Profiling. 
                    PLoS One. 2015;10(9):e0139006. Published 2015 Sep 22. doi:10.1371/journal.pone.0139006
                    This script computes the pivot table (n_query_genes * n_target_genome), 
                    Outputs 1. normalized pivot table (normalized.pivot.npz) 2. discretized pivot table (binned.pivot.npz) to output directory.
                    Query gene id is converted to row number based on (gene_id.pickle); same for target genomes (genome_id.pickle)
                """

    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-r", "--representing", dest="repr",
                  help="FASTA file for representing gene")
    parser.add_option("-d", "--diamond",dest="diamond",
                  help="folder containing diamond blastp results to all reference genome")
    parser.add_option("-b", "--bins",dest="bins", default = 100, type = "int",
                  help="# of bins for discretization")
    parser.add_option("-o", "--outdir",dest="outdir", default = '/tmp', type = "string",
                  help="output directory")

    (options, args) = parser.parse_args()   

    return options

if __name__ == '__main__':
    
    options = option_parser()

    normalized_pivot_table, binned_pivot_table = main(options.repr, options.diamond, bins = options.bins, outdir=options.outdir)

    # save to output directory
    print('results saved in: {}'.format(options.outdir))
    save_npz(os.path.join(options.outdir, 'normalized.pivot.npz'), normalized_pivot_table.tocsr())
    save_npz(os.path.join(options.outdir, 'binned.pivot.npz'), binned_pivot_table)




    