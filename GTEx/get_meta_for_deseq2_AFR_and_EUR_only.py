import pandas as pd
import sys
#upload meta file
metaA = pd.read_csv(sys.argv[1], index_col = 0)

#get an R readable index
metaA.index = [x.replace('-','.') for x in list(metaA.index)]
#define African and European mtDNA haplogroups
afr = list(set([x for x in list(metaA['major_mtDNA_haplogroup']) if x[0] == 'L']))
eur = list(set([x for x in list(metaA['major_mtDNA_haplogroup']) if x[0] != 'L' and x not in ['A2d', 'A2v1', 'B2b3a', 'B4a1a1a18', 'C1c6', 'D4e5b', 'M1a3a', 'M7b1a1a3', 'M7c1c2', 'R32']]))
#assign AFR or EUR to samples based on mtDNA haplorgoup
maternal = []
for i in metaA['major_mtDNA_haplogroup']:
    if i in afr:
        maternal.append('AFR')
    elif i in eur:
        maternal.append('EUR')
    else:
        maternal.append('other')
metaA['maternal'] = maternal
#leave only samples from subjects of African and European origins
metaA = metaA.loc[metaA['maternal'] != 'other']
#save new meta file
metaA.to_csv(sys.argv[2])
