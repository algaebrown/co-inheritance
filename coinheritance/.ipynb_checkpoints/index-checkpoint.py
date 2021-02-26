from os import listdir
from os.path import isfile, join
def gene_id(representing_gene_fasta):
    '''
    replacing representing gene header with shorter index like g0, g1000 to save disk/memory
    :param representing_gene_fasta: .faa containing all representing gene of a pan-genome
    :return: dictionary keys = original gene name; values = id
    '''

    i = 0
    gid = {}
    from Bio import SeqIO
    with open(representing_gene_fasta, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            gid[record.id] = i
            i+=1
    return(gid)

def target_genome_id(target_genome_folder):
    '''
    replacing target_genome_id header with shorter index like t0, t1000 to save disk/memory
    :param target_genome_folder:
    :return: target_id, keys = original genome ID, values = new ID
    '''
    target_lists = [f for f in listdir(target_genome_folder) if isfile(join(target_genome_folder, f))]

    i = 0
    tid = {}
    for target_genome in target_lists:
        tid[target_genome] = i
        i += 1
    return(tid)
