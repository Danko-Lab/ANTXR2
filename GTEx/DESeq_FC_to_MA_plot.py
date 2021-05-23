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
#define the output diractiry
outPutDir = str(sys.argv[2])
#get a list of all files in the defined DESeq2 differential expression by tissue diractiry
fileList = os.listdir(DEdir)

def label_point(x, y, val, ax):
	'''
	Allows for labeling of dots in a scatterplot
	Input:
	x,y - scatterplot coordinates
	val - the label value
	ax - a scatterplot object
	'''
	a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
	for i, point in a.iterrows():
		ax.text(point['x']+0.02, point['y'], str(point['val']), size = 8, weight='bold', style = 'italic', fontname="Arial")

#loop for creating a per-tissue MA plot
for file in fileList:
	print (file)
	plt.figure(file)
	df = pd.read_csv(DEdir + file, names = ['Gene','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']).dropna()
	df[['baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']] = df[['baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']].astype(float)
	color = []
	for p in df['padj'].astype(float):
		if p <= 0.05:
			color.append('red')
		else:
			color.append('gray')
	ax = sns.scatterplot(x = 'baseMean', y = 'log2FoldChange', c = color, alpha = 0.5, data = df, legend = False, edgecolor = 'none', s=10)
	ax.set_xscale('log')
	dfLabel = df.loc[df['Gene'] == 'ENSG00000163297.16']
	dfLabel['Gene'] = ['ANTXR2']
	label_point(dfLabel['baseMean'], dfLabel['log2FoldChange'], dfLabel['Gene'], ax)
	plt.ylabel('Europeans over Africans\nExpression Fold Change (log2)', size = 16)
	plt.xlabel('Mean Expression', size = 16)
	plt.title(file.split('.')[0], size = 40)
	plt.savefig(outPutDir + file.split('.')[0])
