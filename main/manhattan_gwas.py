
import pandas as pd
import numpy as np
from ggplot import *


def gwas_plot(f, 
              column_dict =  {'./gwas_vep.10.tab':'MarkerName', 'p.value':'P','snps':'snps'},
              sorter = np.append(range(1,23,1),['X'])):
    
    data = pd.read_csv(f,sep = '\t',index_col = False,usecols = list(column_dict))
    filename = os.path.basename(f)
    
    data = data.rename(columns = column_dict)
    data = data.dropna(axis=0, how='any')
    data = data.drop_duplicates()
    
    data['CHR'] = data['snps'].str.split('_').str[0]
    data['BP']  = data['snps'].str.split('_').str[1]
    del data['snps']

    print('{} unique snps in {}'.format(data.MarkerName.nunique(), filename))
    
    data['BP'] = data['BP'].astype(int)
    data['P'] = data['P'].astype(float)
    
    chrom_len = data.groupby('CHR')['BP'].max().reset_index() # find max bp on each chrom
    chrom_len['CHR_rank'] = chrom_len['CHR'].map(dict(zip(sorter,range(len(sorter)))))
    chrom_len = chrom_len.sort_values('CHR_rank') # rank the chrom


    BP = chrom_len.BP.tolist()
    cum_BP = [0]
    for i in range(1,len(sorter)):
        cum_BP.append(cum_BP[i-1] + BP[i-1])  
    chrom_len['cum_BP'] = cum_BP
    chrom_len['center'] = chrom_len['cum_BP'] + chrom_len['BP']/2

    cum_data = pd.merge(data,chrom_len, on = 'CHR', how = 'inner',suffixes=('_data', '_max'))
    cum_data['cum_BP'] = cum_data['cum_BP'] + cum_data['BP_data']
    cum_data['-log10P'] = -np.log10(cum_data['P'])
    cum_data["CHR_rank"] = cum_data["CHR_rank"].astype('category')
    cum_data = cum_data[['MarkerName','-log10P','CHR_rank','cum_BP']]


    colors = ["blue", "orange"] * (len(sorter)//2) 
    if len(sorter)%2 != 0:
        colors.append('blue')

    p = ggplot(cum_data,aes(x='cum_BP', y='-log10P', color='CHR_rank')) + \
        geom_point() + \
        scale_color_manual(values = colors) + \
        scale_x_continuous('Chrom',labels = chrom_len['CHR'].tolist(), breaks= chrom_len['center'].tolist()) + \
        ggtitle("{}".format(filename)) 

    p.save('{}.png'.format(filename))
    print('manhathan plot for {}'.format(filename))
    

if __name__ == "__main__":
    gwas_plot(find_file(custom_input['input_dir'],['gwas_sun'])[0])
           
           
    