import pandas as pd
import numpy as np
import os
from ggplot import *

import rpy2 
from rpy2.robjects import r, pandas2ri,numpy2ri
from rpy2.robjects.packages import importr
r['options'](warn=-1)
pandas2ri.activate()
numpy2ri.activate()

pd.set_option('precision', 15)
pd.options.mode.chained_assignment = None
pd.set_option("display.max_rows",100)
pd.set_option("display.max_columns",6000)
pd.set_option('max_colwidth',100)
pd.set_option('display.max_colwidth', -1) 


def chromosome_sort_helper(chrom_list):
    
    chrom_list = [str(chrom) for chrom in chrom_list]
    digits = np.array(sorted([chrom for chrom in chrom_list if chrom.isdigit()]))
    xyz = sorted([chrom for chrom in chrom_list if not chrom.isdigit()])
    
    return np.append(digits, xyz)

def chrom_len_helper(data, CHR = 'CHR'):    
    sorter = chromosome_sort_helper(data[CHR].unique())
    chrom_len = data.groupby(CHR)['BP'].max().reset_index() # find max bp on each chrom
    chrom_len['CHR_rank'] = chrom_len[CHR].map(dict(zip(sorter,range(len(sorter)))))
    chrom_len = chrom_len.sort_values('CHR_rank') # rank the chrom
    
    BP = chrom_len.BP.tolist()
    cum_BP = [0]
    for i in range(1,len(sorter)):
        cum_BP.append(cum_BP[i-1] + BP[i-1])  
    chrom_len['cum_BP'] = cum_BP
    chrom_len['center'] = chrom_len['cum_BP'] + chrom_len['BP']/2
    chrom_len["CHR_rank"] = chrom_len["CHR_rank"].astype('category')
    
    del chrom_len['BP']
    return chrom_len

def eqtl_plot(f, path,
              sep=',',
              column_dict =  {'./eqtl.10.tab':'Gene','variant_id':'snps'}):
    
    filename = os.path.basename(f)
    save_path = os.path.join(path, '{}_eqtl_ggplot'.format(filename))
    if not os.path.isdir(save_path):
        os.makedirs(save_path)
    
    data = pd.read_csv(f,sep = sep,index_col = False,usecols = list(column_dict))
    data = data.rename(columns = column_dict)
    data = data.dropna(axis=0, how='any')
    data['CHR'] = data['snps'].str.split('_').str[0]
    data['BP']  = data['snps'].str.split('_').str[1]
    data['CHR'] = data['CHR'].astype(str)
    data['BP'] = data['BP'].astype(int)
    
    
    print('{} unique snps in {}'.format(data.snps.nunique(), filename))
    print('{} unique Gene in {}'.format(data.Gene.nunique(), filename))
    
    biomaRt = importr('biomaRt') 
    ensembl = biomaRt.useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",GRCh=37)
    gene_data = biomaRt.getBM(
        attributes= np.array(['ensembl_gene_id','chromosome_name','start_position','end_position']),
        filters = 'ensembl_gene_id', values = data.Gene.unique(), mart = ensembl)
    gene_data = pandas2ri.ri2py(gene_data)   
    gene_data['BP'] = (gene_data['start_position'] + gene_data['end_position'])/2
    gene_data['chromosome_name'] = gene_data['chromosome_name'].astype(str)
    
    gene_chrom_len = chrom_len_helper(gene_data,CHR = 'chromosome_name')
    gene_data_1 = pd.merge(gene_data,gene_chrom_len, on = 'chromosome_name', how = 'inner')
    gene_data_1['cum_BP'] = gene_data_1['cum_BP'] + gene_data_1['BP']
    gene_data_1 = gene_data_1[['ensembl_gene_id','cum_BP']]

    snps_chrom_len = chrom_len_helper(data)
    snps_data = pd.merge(data,snps_chrom_len, on = 'CHR', how = 'inner')
    snps_data['cum_BP'] = snps_data['cum_BP'] + snps_data['BP']
    snps_data = snps_data[['Gene','snps','cum_BP']]
    
    plot_data = pd.merge(snps_data, gene_data_1,
                     left_on = 'Gene',
                     right_on = 'ensembl_gene_id',
                     how = 'inner',
                     suffixes = ('_snps','_gene'))
    gene_chrom_len.to_csv(os.path.join(save_path,'gene_chrom_len.csv'), index = False)
    snps_chrom_len.to_csv(os.path.join(save_path,'snps_chrom_len.csv'), index = False)
    plot_data.to_csv(os.path.join(save_path,'plot_data.csv'), index = False)

    p = ggplot(plot_data,aes(x='cum_BP_snps', y='cum_BP_gene')) + \
        geom_point() + \
        scale_x_continuous('snps',labels = snps_chrom_len['CHR'].tolist(), breaks= snps_chrom_len['center'].tolist()) + \
        scale_y_continuous('gene',labels = gene_chrom_len['chromosome_name'].tolist(), breaks= gene_chrom_len['center'].tolist()) + \
        ggtitle("{}".format(filename)) 

    p.save(os.path.join(save_path,'{}.png'.format(filename)))
    print('manhathan plot for {}'.format(filename))

if __name__ == "__main__":
    eqtl_plot(eqtl_file, save_path, column_dict = {'gene_id':'Gene','variant_id':'snps'})