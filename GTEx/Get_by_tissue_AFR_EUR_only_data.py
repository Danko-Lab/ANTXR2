import pandas as pd
import sys
# get AFR and EUR samples raw reads
reads = pd.read_csv(sys.argv[1], sep = '\t', index_col = 0)
# read in the sample-tissue assingment from GTEx
sampleInfo = pd.read_csv(sys.argv[2], usecols = [0,2], names = ['SAMPID', 'SMTSD'])
outputDir = str(sys.argv[3])
sampleInfo['SAMPID'] = sampleInfo['SAMPID'].astype(str)
#make sure R-adjusted raw reads sample names match the ones on the GTEx sample-tissue assingment file
sampleInfo['SAMPID'] = [x.replace('-','.') for x in list(sampleInfo['SAMPID'])]
#make sure the new raw reads file only contains samples that are in the GTEx sample-tissue assingment file
samples_to_keep = [x for x in list(reads.columns) if x in list(sampleInfo['SAMPID'])]
lmData = reads[samples_to_keep].transpose()
#assign tissue of sample origin to the raw reads file's samples
tissue = []
c = 0
for sample in samples_to_keep:
	sampleInfo_sample = sampleInfo.loc[sampleInfo['SAMPID'] == sample]
	tissue.append(list(sampleInfo_sample['SMTSD'])[0])
	c += 1
	print ('Assigning sample ' + str(c) + ' out of ' + str(len(samples_to_keep)))
lmData['tissue'] = tissue
#split the raw reads file by tissue to get a per-tisse files directory
for tis in list(set(tissue)):
	lmData.loc[lmData['tissue'] == tis][lmData.columns[:-1]].transpose().to_csv(sys.argv[3] + str(tis) + '_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads_AFR_EUR_only_for_DEseq2.gct.gz', sep = '\t', compression = 'gzip')



