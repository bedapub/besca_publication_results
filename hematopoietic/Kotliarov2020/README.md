# Kotliarov Dataset

The publicly available data was processed using the extract_data.R script to extract required data to convert into mtx format for loading into besca.  
This data was then loaded and parsed into the correct besca format using parsing_public_data.ipynb.  
The raw data was processed using the besca standard workflow as outlined in standard_workflow_besca_2.0.ipynb.  
The clr_normalized adata object generated by the besca standard workflow was used as input to the R script convert_scRNAseq_to_fcs.R which generated an .fcs file which could be loaded into FlowJo and manually gated there.  
The gated cell popualtions were exported to .fcs files from FlowJo and converted back to lists which can be loaded into scanpy using the R script FlowJo_extract_barcodes.R. 
Celltype annotation was performed in the jupyter notebook celltype_annotation_besca.ipynb.  

