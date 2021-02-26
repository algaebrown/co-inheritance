import os
from Bio import SeqIO
import pickle

def gene_id(representing_gene_fasta, outdir = '.'):
    '''
    replacing representing gene header with shorter index like g0, g1000 to save disk/memory
    :param representing_gene_fasta: .faa containing all representing gene of a pan-genome
    :return: dictionary keys = original gene name; values = id
    '''
    full_path = os.path.join(outdir, 'gene_id.pickle')
    
    if os.path.isfile(full_path):
        print('Exist gene index {}'.format(full_path))
        with open(full_path, 'rb') as handle:
            gid = pickle.load(handle)
    else:
        print('Writing gene index file to {}'.format(full_path))
        i = 0
        gid = {}
    
        with open(representing_gene_fasta, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                gid[record.id] = i
                i+=1
        with open(full_path, 'wb') as handle:
            pickle.dump(gid, handle)
    
    return(gid)

def target_genome_id(target_genome_folder, outdir = '.'):
    '''
    replacing target_genome_id header with shorter index like t0, t1000 to save disk/memory
    :param target_genome_folder:
    :return: target_id, keys = original genome ID, values = new ID
    '''
    
    full_path = os.path.join(outdir, 'genome_id.pickle')
    
    if os.path.isfile(full_path):
        print('Exists genome index {}'.format(full_path))
        with open(full_path, 'rb') as handle:
            tid = pickle.load(handle)
    else:
        print('Writing genome index file to {}'.format(full_path))
        target_lists = [f for f in os.listdir(target_genome_folder) if os.path.isfile(os.path.join(target_genome_folder, f))]
    
        i = 0
        tid = {}
        for target_genome in target_lists:
            tid[target_genome] = i
            i += 1
        with open(full_path, 'wb') as handle:
            pickle.dump(tid, handle)
    return(tid)
