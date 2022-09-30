## runafter merge_gene

import os
import sys
import pandas as pd
import numpy as np
import warnings
import shutil
import collections
import pickle


from filef import find_file,create_folder
from generate_custom_input import open_custom_input, save_custom_input



def read_gene(files):
    result = pd.DataFrame()
    for f in files:
        df_csv = pd.read_csv(f,index_col = False)
        df_csv['partition'] = f.split('/')[-1]
        result = pd.concat([result,df_csv])
    return result

def original_t_analysis(data):

    z_eqtl = np.abs(data.Z_eqtl) 
    z_gwas = np.abs(data.Z_gwas)
    min_comp = list(np.minimum(z_eqtl, z_gwas)) # min z of each SNP
    t = max(min_comp) if min_comp else None   # t value as max of z list
    
    t_table = collections.defaultdict(list)
    if t:
        t_table['t'].append(t)

        # count how many z in eqtl file bigger than t
        m = sum(z_eqtl >= t)  
        t_table['m_eqtl'].append(m) 

        # count how many z in summary file bigger than t
        k = sum(z_gwas >= t)
        t_table['k_gwas'].append(k) 

        index = min_comp.index(t)
        if  z_eqtl.iloc[index] == t:
            t_table['t_from'].append('eqtl')
        else:
            t_table['t_from'].append('gwas')

        fill = data.iloc[index]

        t_table['MarkerName_given_t'].append(fill['MarkerName'])
        
        t_table['Z_eqtl_for_MarkerName_given_t'].append(fill['Z_eqtl'])
        t_table['Beta_eqtl_for_MarkerName_given_t'].append(fill['Beta_eqtl'])
        t_table['Z_gwas_for_MarkerName_given_t'].append(fill['Z_gwas'])
        t_table['Beta_gwas_for_MarkerName_given_t'].append(fill['Beta_gwas']) 
        
        t_table['Gene'].append(fill['Gene'])
        t_table['Snps_given_t'].append(fill['Snps'])
        
        t_table = pd.DataFrame(t_table)
        return t_table


if __name__ == "__main__":
    
    result_name = sys.argv[1]
    node = int(sys.argv[2])
    print('work on node {}'.format(node))
    
    cwd = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))
    
    result_dir = os.path.join(cwd,result_name) 
    gene_dir = create_folder(result_dir,'preprocessed_gene')
    gene_dir_2 = create_folder(result_dir,'gene_merged')
    
    t_dir1 = create_folder(result_dir,'t_org')
      
    geneD = open_custom_input(os.path.join(result_dir,'geneD.txt'))
    parralle_geneD = open_custom_input(os.path.join(result_dir,'parralle_geneD.txt'))    
    print('{} genes in this node'.format(len(parralle_geneD[node])))
    
    
    result1 = pd.DataFrame()
    for gene in parralle_geneD[node]:
        
        # if not os.path.isfile(os.path.join(gene_dir_2,'{}.csv'.format(gene))):
        data1 = read_gene(geneD[gene])
        data1 = data1.drop_duplicates()
        print('{} snps in {}'.format(data1.Snps.nunique(), gene))
        
        data1.to_csv(os.path.join(gene_dir_2,'{}.csv'.format(gene)),index = False)
        result1 = pd.concat([result1,original_t_analysis(data1)])
        print('done merge {}'.format(gene))
        
    
    result1.to_csv(os.path.join(t_dir1,'t_{}.csv'.format(node)),index = False) 
    print(os.path.join(t_dir1,'t_{}.csv'.format(node)))
    
    
    
