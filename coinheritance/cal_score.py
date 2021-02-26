from scipy.sparse import load_npz
from itertools import combinations
from sklearn.metrics import mutual_info_score, normalized_mutual_info_score
import numpy as np
from multiprocessing import Process, Manager
import multiprocessing
import time
import pandas as pd
import os


def get_mutual_info(index_1, index_2,  score = 'mutual_info', unique_index = False):
    ''' given two row index(genes), extract information from binned matrix and calculate mutual informations
    specify score='normalized_mutual_info' for normalized mutual information
    returns int
    '''
    row_1 = binned[index_1, :].toarray().ravel()
    row_2 = binned[index_2, :].toarray().ravel()
    
    if np.sum(row_1) == 0 or np.sum(row_1) == 0:
        # no data
        return 0
    
    if score == 'mutual_info':
        s = mutual_info_score(row_1, row_2)
    elif score == 'normalized_mutual_info':
        s = normalized_mutual_info_score(row_1, row_2)
    if not unique_index:
        results.append((index_1, index_2, s)) # save to common memory
    else:
        # make sure small index goes first
        if index_1<index_2:
            results.append((index_1, index_2, s))
        else:
            results.append((index_2, index_1, s))
    return s

def option_parser():
    from optparse import OptionParser
    usage = """
        THIS IS CO-INHERITANCE 1.0.0
        python cal_score.py --table ~/data0118/binned.pivot.npz -o ~/data0118
        """
    description = """CO-INHERTITANCE is a tool to compute co-inheritance relationship for a set of genes.
                    Methods is by this paper:
                    Shin J, Lee I. Co-Inheritance Analysis within the Domains of Life Substantially Improves Network Inference by Phylogenetic Profiling. 
                    PLoS One. 2015;10(9):e0139006. Published 2015 Sep 22. doi:10.1371/journal.pone.0139006
                    This script computes the mutual information/normalized mutual information based on pivot table (output by diamond_to_pivot.py)
                """

    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-p", "--pool", dest="pool", type = "int", default = 8,
                  help="number of multiprocessing threads")
    parser.add_option("-t", "--table",dest="pivot",
                  help="path to pivot table (binned.pivot.npz)")
    parser.add_option("-s", "--score",dest="score", default = 'mutual_info', type = "string",
                  help="mutual_info or normalized_mutual_info")
    parser.add_option("-o", "--outdir",dest="outdir", default = '/tmp', type = "string",
                  help="output directory")
    parser.add_option("-n", "--subset",dest="n_gene", default = 0, type = "int",
                  help="select n genes to do computation; debugging use or test run")

    (options, args) = parser.parse_args()   

    return options

if __name__=='__main__':
    options = option_parser()

    binned = load_npz(options.pivot)

    if options.n_gene > 0:
        # debugging option
        binned = binned[:options.n_gene, :]

    t1 = time.time()
    # only run on genes with hits in the reference genome
    non_zero_genes = np.where(np.sum(binned, axis = 1)>0)[0] 
    all_index_combine = combinations(list(non_zero_genes), 2)
    tasks = list(all_index_combine)
    
    print('No. tasks: {}'.format(len(tasks)))
    # shared memory for saving results
    manager = Manager()
    results = manager.list()

    for i in range(len(tasks)):
        keywords = {'score':options.score}
        p = Process(target=get_mutual_info, args=tasks[i], kwargs=keywords)
        p.start()
        p.join()
    
    t2 = time.time()
    print('Done computing scores, extracting results to dataframe. Total time: {} secs'.format(t2-t1))
    
    df = pd.DataFrame(list(results), columns = ['gene_one', 'gene_two', 'mutual_info'])

    df.to_csv(os.path.join(options.outdir, options.score+'.csv'), index = False)

    print('results saved to: {}'.format(os.path.join(options.outdir, options.score+'.csv')))
