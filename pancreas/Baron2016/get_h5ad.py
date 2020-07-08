## THIS SCRIPT WAS USED TO CREATE A CONCATENATED RAW H5AD files starting with individual samples download from GEO.
import scanpy as sc
import besca as bc
import pandas as pd

folder_to_sample = ['GSM2230757_human1_umifm_counts',  'GSM2230758_human2_umifm_counts' , 'GSM2230759_human3_umifm_counts', 'GSM2230760_human4_umifm_counts']

sample_name = [x.split('_')[1] for x in folder_to_sample]
sample_name

import pandas as pd
def read_in( name ):
    indir='downloaded/' + name
    print(indir)
    # GET FULL MATRIX
    newFolder = 'raw_samples/' + name
    counts = pd.read_csv(indir +'.csv',header=0,sep=',', index_col=0)
    # GET METADATA
    metadata = counts.loc[:,['barcode','assigned_cluster']]
    # REMOVE METADATA
    counts = counts.drop(['barcode','assigned_cluster'], axis=1).transpose()
   
    return metadata, counts



import os
def newFF(newFolder):
    if not os.path.exists(newFolder):
        os.mkdir(newFolder)
        print("Directory " , newFolder ,  " Created ")
    else:    
        print("Directory " , newFolder ,  " already exists")




for name in folder_to_sample:
    newFolder = 'raw_samples2/' + name
    newFF(newFolder)
    metadata, counts  = read_in( name )    
    metadata.head()
    genes  = counts.reset_index()['index']
    genes.to_csv(newFolder+ '/genes.tsv', sep = '\t', index = False)
    counts.transpose().to_csv(newFolder + '/counts.tsv',sep='\t')
    metadata.to_csv(newFolder + '/metadata.tsv')# TO WRITE
    newFolder = 'raw_samples2/' + name
    aa = sc.read_csv(newFolder + '/counts.tsv',delimiter='\t')
    metadata = pd.read_csv(newFolder + '/metadata.tsv', index_col= 0)
    tmp = pd.merge( aa.obs , metadata, how='left', left_index=True, right_index=True)
    aa.obs  = tmp
    aa.write_h5ad('raw/' + name + '.h5ad')



# NAME
adatas = []
for name in folder_to_sample:
    print(name)
    adatas += [sc.read_h5ad('raw/' + name + '.h5ad') ]


final_adata = adatas[0].concatenate(*adatas[1:],   
                                     batch_key= 'Individual', batch_categories = folder_to_sample)


final_adata.write_h5ad('raw/all_samples.h5ad')


####### NEXT

import scanpy as sc
aa = sc.read_h5ad( 'raw/all_samples.h5ad')
aa.var.to_csv('raw/genes_Symbol.csv')