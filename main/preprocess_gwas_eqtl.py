##it is an illustration example in SSSD_2022b for dataset SummaryStatistics.csv and eqtl from AD 2016 

import os
import sys
import pandas as pd
import numpy as np
import warnings

import collections
import pickle


from filef import find_file,create_folder
from generate_custom_input import open_custom_input, save_custom_input

def create_eqtl(file):
    tb = pd.read_table(file, index_col = False)
    tb.columns = ['Snps','Gene','Z_eqtl','Pvalue_eqtl',
                  'FDR_eqtl','Beta_eqtl']
    return tb

def merge_eqtl(dir_name, chron):

    cis = os.path.join(dir_name, 'eqtl.' +  str(chron) + '.cis.tab')
    trans = os.path.join(dir_name, 'eqtl.' + str(chron) + '.trans.tab')
    
    t1 = create_eqtl(cis)
    t2 = create_eqtl(trans)       
    
    eqtl = pd.concat([t1, t2],ignore_index = True)
    eqtl = eqtl.drop_duplicates(subset = ['Snps','Gene'])
    return eqtl

def create_gwas(file):
    tb = pd.read_table(file, index_col = False,sep = ',')

    tb.columns = ['MarkerName',
         'Pvalue_gwas',
         'chromosome_build38',
         'base_pair_location_build38',
         'CHR',
         'POS',
         'Effect_allele',
         'Non_Effect_allele',
         'effect_allele_frequency',
         'odds_ratio',
         'ci_lower',
         'ci_upper',
         'Beta_gwas',
         'SE_gwas',
         'n_cases',
         'n_controls',
         'het_isq',
         'het_pvalue',
         'variant_alternate_id']
    
    tb['Beta_gwas'] = pd.to_numeric(tb['Beta_gwas'],errors='coerce')
    tb['SE_gwas']   = pd.to_numeric(tb['SE_gwas'],errors='coerce')
    
    tb['Z_gwas'] = tb.Beta_gwas.astype('float64')/tb.SE_gwas.astype('float64')
    tb = tb.dropna(axis=0, how='any')
    return tb



if __name__ == "__main__":
        
    cwd = '/home/shuopwu/run/SSSD_2022'
    # cwd = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) #project dir

    ##put raw data into the data dir 
    data_dir = os.path.join(cwd,'data')
    eqtl_dir = os.path.join(data_dir,'SSSD_AD_2016')
    gwas_file = os.path.join(data_dir,'SummaryStatistics.csv')

    result_dir = os.path.join(cwd,'result') 
    preprocessed_dir = create_folder(result_dir,'preprocessed_data')
    
    ##preprocessing data for desired column names and snps subset 
    gwas_data = create_gwas(gwas_file)

    for i in range(22):
        chron = i + 1
        e = merge_eqtl(eqtl_dir,chron)
        g = gwas_data[gwas_data['CHR'] == chron]
        
        ##reconstruct snps as chr:pos_a1_a2
        g['Snps'] = g[['CHR', 'POS']].astype(str).agg(':'.join, axis=1) # handle NaNs
        g['Non_Effect_allele'] = g['Non_Effect_allele'].str.upper()
        g['Effect_allele'] = g['Effect_allele'].str.upper()
        g['Snps'] = g[['Snps','Non_Effect_allele','Effect_allele']].agg('_'.join, axis=1)  
            
        data = pd.merge(left = e, right = g, how = 'inner', on = 'Snps',sort = False)
        
        
        print('-----------------------------------------------')
        print('for chr_{}'.format(chron))
        print('in eqtl {} Snps and {} Gene'.format(e.Snps.nunique(),e.Gene.nunique())) 
        print('in gwas {} Snps'.format(g.Snps.nunique()))
        print('we merged {} Snps and {} Gene'.format(data.Snps.nunique(),data.Gene.nunique()))
        print('-----------------------------------------------')

        data.to_csv(os.path.join(preprocessed_dir,'chr_{}.csv'.format(chron)),index = False)
        print(os.path.join(preprocessed_dir,'chr_{}.csv'.format(chron)))
    
    

   