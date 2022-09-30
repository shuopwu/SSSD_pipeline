########################
### basic functions ####
########################

import os
import sys
import pandas as pd
import numpy as np
import warnings
import collections
import shutil
# from pathlib import Path
from datetime import datetime

from filef import find_file,create_folder
from generate_custom_input import open_custom_input, save_custom_input


##below is for bsub<.sh 
# def hpc_load(sh,name,error_dir):
#     with open(sh, 'w') as f1: 
#         f1.write('''
# #!/bin/bash
# #BSUB -J {}            # LSF job name
# #BSUB -o {}/{}.out     # Name of the job output file 
# #BSUB -e {}/{}.error   # Name of the job error file
# #BSUB -M 20480

# source $HOME/env/my_python-3.6.3/bin/activate
# module load R/4.0.2

# '''.format(name,error_dir,name,error_dir,name))
#     print('{} generated'.format(sh))    

##below is for the bsub per line
def hpc_line(name,error_dir,wait = [],memory = 10240):
    temp = 'bsub -J {} -o {}/{}.out -e {}/{}.error -M {} '.format(name,error_dir,name,error_dir,name,memory)
    if len(wait) != 0:
        temp = temp + '-w "'
        for i in wait:
            temp = temp + 'done({}) && '.format(i)
        temp = temp[:-4] + '" '
    return temp

if __name__ == "__main__":
    ##eg. python bash_generate.py MFG 20 
    disease = sys.argv[1]
    nodes = int(sys.argv[2]) # used for parallel, usually try 20 
    
    cwd = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) #project dir 
    result_name = 'result_{}_{}'.format(disease,datetime.now().strftime('%m%d'))
    
    print('project path is {}'.format(cwd))
    
    
    ###################################################
    ############ automatic generated below ############
    ###################################################
    
    '''
    # below is when inside the interactivenode envs 
    # bsub_M # alias = bsub -Is bash
    # python this file
    # sh $script
    '''
    
    ##paralle genes preprocessing
    # hpc_load(script,'run','error') #no need to do this if already inside the interactive node envs
       
    ##functions files 
    function_dir = create_folder(cwd,'main')
    # bed_dir = os.path.join(cwd,'jianqiao')
    
    ##result dir 
    result_dir = create_folder(cwd,result_name)
    # dirpath = Path(result_dir) / 'error'
    # if dirpath.exists() and dirpath.is_dir():
    #     shutil.rmtree(dirpath)
    error_dir = create_folder(result_dir,'error')
    tmp_dir = create_folder(result_dir,'tmp')
    
    
    ##script file 
    script = os.path.join(result_dir,'run.sh')
    # script_dir = create_folder(result_dir,'node_scripts')
    
    generate_gene_f = os.path.join(function_dir,'generate_gene.py')
    gene_dir = create_folder(result_dir,'preprocessed_gene')
    
    merge_gene_f = os.path.join(function_dir,'merge_gene.py')
    t_f = os.path.join(function_dir,'t.py')
 
 
    ############################################################ 
    ############ manually put the project dir below ############
    ############################################################ 
    data_dir = create_folder(cwd,'data') #if you want to upload later into the project cwd
    # data_dir =  '/home/shuopwu/data/SSSD_2022' #if data dir is not inside the project cwd
     
    ##upload preprocessed data to result folder if already preprocessed in local,and was uploaded from local to data dir 
    if os.path.isdir(os.path.join(data_dir,'preprocessed_data')) and not os.path.isdir(os.path.join(result_dir,'preprocessed_data')):
                       
        shutil.copytree(os.path.join(data_dir,'preprocessed_data'),
                        os.path.join(result_dir,'preprocessed_data'))
        shutil.rmtree(os.path.join(data_dir,'preprocessed_data'))
        
    preprocessed_dir = create_folder(result_dir,'preprocessed_data')
        
    
    # hpc_load(script,'run','error') #no need if already in bsub_I envs 
    with open(script , 'w') as f5:
        
        f5.write('#!/bin/bash')
        f5.write('\n')
        f5.write('cd {}'.format(cwd)) #cwd is the project dir 
        f5.write('\n')
        
        ###################################
        ###### preprocessed_raw_data ######
        ###################################
        
        ##preprocessing in pipleline, usually we run this in local envs and upload the preprocessed_data to  preprocessed_dir
        ##see Notion (SSSD 2016 Processing Step by Step review) as an example 
        ##see preprocessing_MFG_chunks.py in Notion (SSSD 2022 Processing Step by Step review) as an example 
        ##see merge_gwas_eqtl.py in this folder as an example 
        
        # f5.write(hpc_line('{}_chunk'.format(disease),'error',
        #                       memory = 20480) + 'python main/preprocessing_MFG_chunks.py {}'.format(disease))
        # f5.write('\n')
        print('preprocessed_gwas_eqtl loaded')  
        
        ###################################
        ###### preprocessed_gene ##########
        ###################################
        wait_list_preprocessed_gene = []
        
        preprocessed_files = find_file(preprocessed_dir)
        for i in range(len(preprocessed_files)):
            f = preprocessed_files[i]
            chunk_name = os.path.basename(f).split('_')[1][:-4]#.split('.')[0]

            ##pipeline with wait for preprocessing raw data, usually we do preprocessing raw data in local and skip this part 
            # f5.write(hpc_line('generate_gene_{}_{}'.format(chrom_chunks,disease),'error',
            #                   wait = ['{}_chunk'.format(disease)],
            #                   memory = 20480) + 'python {} {} {}\n'.format(generate_gene_f,result_name,chrom_chunks))
            
            
            ## single module 
            f5.write(hpc_line('generate_gene_{}_{}'.format(chunk_name,disease),'error',memory = 20480) + 'python {} {} {}\n'.format(
                                  generate_gene_f,result_name,chunk_name))
            
            wait_list_preprocessed_gene.append('generate_gene_{}_{}'.format(chunk_name,disease))
            
            
        print('preprocessed_gene loaded')  
        f5.write('sleep 5s')
        f5.write('\n')
        
        ###################################
        ###### merge_gene #################
        ###################################
        
        ## pipeline with wait for above step 
        f5.write(hpc_line('merge_gene_{}'.format(disease),'error',wait = wait_list_preprocessed_gene,memory = 20480) + 'python {} {} {}\n'.format(merge_gene_f,result_name,nodes))
        
        ## single module 
        # f5.write(hpc_line('merge_gene','error',memory = 20480) + 'python {} {} {}\n'.format(merge_gene_f,result_name,nodes))
        
        print('merge_gene loaded')  
        f5.write('sleep 5s')
        f5.write('\n')

        ##############################
        ####### t analysis ###########
        ##############################
        for node in range(nodes):
            ## pipeline with wait for above step 
            f5.write(hpc_line('t_{}_{}'.format(node,disease),'error',wait = ['merge_gene_{}'.format(disease)],memory = 20480) + 'python {} {} {}\n'.format(t_f,result_name,node))
             
            ## single module      
            # f5.write(hpc_line('t_{}'.format(node),'error',memory = 20480) + 'python {} {} {}\n'.format(t_f,result_name,node))
            
        print('t_analysis loaded')  
        f5.write('sleep 5s')
        f5.write('\n')
        
        #########################################
        ####### Jianqiao Sat analysis ###########
        #########################################
        # sat_f1 = os.path.join(function_dir,'sat.py')
        # sat_f2 = os.path.join(function_dir,'sat_hpc.R')

        # for node in range(nodes):
        #     ##run in pipleline
        #     f5.write(hpc_line('sat_py_{}_{}'.format(node,disease),'error',wait = ['t_{}_{}'.format(node,disease)],memory = 20480) + 'python {} {} {}\n'.format(sat_f1,result_name,node))

            ##run below in R only after you have test run jianqiao's code successfully 
            # f5.write(hpc_line('sat_r_{}_{}'.format(node,disease),'error',wait = ['sat_py_{}_{}'.format(node,disease)],memory = 20480) + 'Rscript {} {} {} {}\n'.format(sat_f2,result_name,node,cwd))           
            
        # f5.write('\n')
        # print('sat_analysis loaded')  
        # f5.write('\n')
        
        print('sh {} '.format(script))
    


    
    