# for hpc only 
# copy from New_Data_Preprocessing_071222
# for dataset SummaryStatistics.csv only 

import os
import sys
import pandas as pd
import numpy as np
import warnings


from filef import find_file,create_folder

if __name__ == "__main__":

    result_name = sys.argv[1]
    chunk_name = int(sys.argv[2]) 
        
    cwd = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))

    result_dir = os.path.join(cwd,result_name) 
    preprocessed_dir = create_folder(result_dir,'preprocessed_data')
    gene_dir = create_folder(result_dir,'preprocessed_gene')

    data = pd.read_csv(os.path.join(preprocessed_dir,'chr_{}.csv'.format(chunk_name)),index_col = False)
    for gene in data.Gene.unique():
        sub_data = data[data.Gene == gene]
        sub_data.to_csv(os.path.join(gene_dir,'{}_chr_{}.csv'.format(gene,chunk_name)),index = False)
        print(gene)
    
    

   