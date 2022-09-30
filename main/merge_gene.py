# merge gene_chr1, gene_chr2,.. gene_chr22 to one file for the same gene
import os
import sys
import pandas as pd
import numpy as np
import warnings

import collections
import pickle


from filef import find_file,create_folder
from generate_custom_input import open_custom_input, save_custom_input



def create_geneD(files):
    geneD = collections.defaultdict(list)
    for f in files:
        gene = os.path.basename(f).split('_')[0]
        geneD[gene].append(f)
        # if '.csv' in f:
        #     
        # if "ENSG" in gene :
        #     

    return geneD

def create_parralle_geneD(geneD, chunks = 10):
    parralle_geneD = collections.defaultdict(dict)
    
    gene_list = np.asarray(sorted(list(geneD.keys())))
    
    if len(gene_list) > chunks:
        slicing = gene_list.shape[0] % (chunks - 1)
        parralle_geneD_list = [gene_list[:slicing]]+ np.split(gene_list[slicing:],chunks - 1)
    
    else: 
        parralle_geneD_list= [gene_list]
                          
    parralle_geneD_list = [gene_list for gene_list in parralle_geneD_list if len(gene_list)!=0]
    
    for i in range(len(parralle_geneD_list)):
        parralle_geneD[i] = parralle_geneD_list[i]
        
    return parralle_geneD

if __name__ == "__main__":
    
    result_name = sys.argv[1]
    nodes = int(sys.argv[2])
    
    cwd = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))
    
    result_dir = os.path.join(cwd,result_name) 
    gene_dir = create_folder(result_dir,'preprocessed_gene')
    
    ##parallel genes merge
    gene_files = find_file(gene_dir)    
    geneD = create_geneD(gene_files)
    parralle_geneD = create_parralle_geneD(geneD, chunks = nodes)
    
    save_custom_input(geneD,os.path.join(result_dir,'geneD.txt'))
    save_custom_input(parralle_geneD,os.path.join(result_dir,'parralle_geneD.txt'))
    
    
    print('done write parralle_geneD, current nodes = {}'.format(nodes))
    print(os.path.join(result_dir,'geneD.txt'))
    print(os.path.join(result_dir,'parralle_geneD.txt'))

    
          
 
