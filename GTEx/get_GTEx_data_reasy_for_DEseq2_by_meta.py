import pandas as pd
import sys
#upload GTEx raw reads count per gene
reads = pd.read_csv(sys.argv[1], sep = '\t', skiprows = [0,1], index_col = 0)
#remove the "Description" column
reads = reads.drop('Description', axis=1)
#change "-" to "." in samples names, to prevent R errors when running DESeq2
reads.columns = [x.replace('-','.') for x in list(reads.columns)]
#read meta file
meta = pd.read_csv(sys.argv[2], index_col = 0)
#make sure new meta will only contain samples in the raw count file
#this is unlikely to result in any loss of samples
meta = meta.loc[meta.index.isin(list(reads.columns))]
#make sure new raw reads file will only contain samples not in the meta file
#this results in a loss of some samples to which we didn't have mtDNA haplogroup data for
reads = reads[list(meta.index)]
meta = meta.transpose()[list(reads.columns)].transpose()
#save new raw reads and meta files
reads.to_csv(sys.argv[3], sep = '\t', compression = 'gzip')
meta.to_csv(sys.argv[4])