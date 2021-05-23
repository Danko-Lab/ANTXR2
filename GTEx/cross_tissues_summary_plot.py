import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as plt
import pandas as pd
import os
import seaborn as sns
import scipy.stats as stats
import numpy as np
import sys
#define the DESeq2 differential expression by tissue diractiry
DEdir = str(sys.argv[1])
#get a list of all files in the defined DESeq2 differential expression by tissue diractiry
fileList = os.listdir(DEdir)
#get a list of all tissues based on file names in the DESeq2 differential expression by tissue diractiry
tissueList = [x.split('.')[0] for x in fileList]
#create a new dataframe to be fed with the DE values and statistics for ANTXR2 (ENSG00000163297.16)
#and loop over the tissues to feed it
data = pd.DataFrame(columns = ['baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'])
c = 0
for file in fileList:
	print (file)
	#plt.figure(file)
	df = pd.read_csv(DEdir + file, names = ['Gene','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']).dropna()
	df[['baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']] = df[['baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']].astype(float)
	antxr2 = df.loc[df['Gene'] == 'ENSG00000163297.16']
	values = [list(antxr2['baseMean'])[0], list(antxr2['log2FoldChange'])[0], list(antxr2['lfcSE'])[0], list(antxr2['stat'])[0], list(antxr2['pvalue'])[0], list(antxr2['padj'])[0]]
	print (values)
	data.loc[c] = values
	c += 1
print (data)
data = data.dropna()
#add the tisse list as a column in the new dataframe
data['Tissue'] = tissueList
#save new per-tissue ANTXR2 specific DE statistics dataframe to csv
data.to_csv(sys.argv[2])
