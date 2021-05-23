import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#get the list of ANTXR2 eSNPs form the DICE dataset
esnps = list(set(['rs28830602', 'rs28830602', 'rs28842234', 'rs28842234', 'rs4405963', 'rs4405963', 'rs66968950', 'rs66968950', 'rs28763597', 'rs28763597', 'rs12499307', 'rs12499307', 'rs35263215', 'rs35263215', 'rs4389526', 'rs4389526', 'rs4640621', 'rs4640621', 'rs4637335', 'rs4637335', 'rs7687189', 'rs7687189', 'rs4333130', 'rs4333130', 'rs4355335', 'rs4355335', 'rs59592559', 'rs72655081', 'rs17510948', 'rs17464552', 'rs72655084', 'rs79168761', 'rs143385483', 'rs79682523', 'rs72655089', 'rs72655092', 'rs72655094', 'rs56984496', 'rs72655096', 'rs61305264', 'rs72655097', 'rs115397861', 'rs72656506', 'rs111402170', 'rs72656507', 'rs114184382', 'rs72656511', 'rs76369210', 'rs72656513', 'rs72656517', 'rs72656521', 'rs3834215', 'rs72656525', 'rs146379067', 'rs116797028', 'rs17465218', 'rs72656528', 'rs72656529', 'rs1458044', 'rs41279283', 'rs72656531', 'rs56102626', 'rs142234659', 'rs13110335', 'rs72656538', 'rs72656541', 'rs372041304', 'rs114563177', 'rs369055754', 'rs78545693']))
#loop over the PLINK-calculated LD scores between the ANTXR2 eSNPs
#and make a directory with the calculated R-sq and D'
columns = []
drdict = {}
filepath = os.path.join('C:\Users','barsh','Desktop','Gilad','ANTXR2','pairwise_LD_results_D_R-sq.txt')
with open(filepath) as fp:
    c = 0
    for line in fp:
        c += 1
        if line[0:2] == 'rs':
            i = 0
            col = ''
            while line[i] != ',':
                col = col + line[i]
                i += 1
            ind = line[i + 2: -1]
            if col not in drdict:
                drdict[col] = {}
                columns.append(col)
            if col in drdict and ind not in drdict[col]:
                nx = next(fp)
                drdict[col][ind] = {'R-sq': float(nx[7:12]), 'D': float(nx[19:-1])}
        c += 1
        print (c)
#createa all 1 dataframes to feed D and R-sq into between all pairs of ANTXR2 eSNPs
df_Rsq_to_heatmap = pd.DataFrame(1.0, columns = columns, index = columns)
df_D_to_heatmap = pd.DataFrame(1.0, columns = columns, index = columns) 
#feed R-sq and D' values
for co in df_Rsq_to_heatmap.index:
    for ix in drdict[co]:
        if df_Rsq_to_heatmap[co][ix] == 1.0 and df_Rsq_to_heatmap[ix][co] == 1.0:
            df_Rsq_to_heatmap[co][ix]  = float(drdict[co][ix]['R-sq'])
            df_Rsq_to_heatmap[ix][co]  = float(drdict[co][ix]['R-sq'])
            df_D_to_heatmap[co][ix]  = float(drdict[co][ix]['D'])
            df_D_to_heatmap[ix][co]  = float(drdict[co][ix]['D'])
        
#draw clustered heatmap based on R-sq values
f, ax = plt.subplots(figsize=(20, 20))
sns.heatmap(df_Rsq_to_heatmap, cmap = sns.color_palette("Reds"), ax=ax, square = True, linewidth = 0)
sns.set(font_scale=1)
f.savefig(os.path.join('C:\Users','barsh','Desktop','Gilad','ANTXR2','pairwise_LD_results_Rsq_heatmap.svg'))

#draw clustered heatmap based on D' values
f, ax = plt.subplots(figsize=(20, 20))
sns.heatmap(df_D_to_heatmap, cmap = sns.color_palette("Reds"), ax=ax, square = True, linewidth = 0)
sns.set(font_scale=1)
f.savefig(os.path.join('C:\Users','barsh','Desktop','Gilad','ANTXR2','pairwise_LD_results_D_heatmap.svg'))


#get the DICE data for the ANTXR2 eSNPs
esnp_data = pd.read_csv(os.path.join('C:\Users','barsh','Desktop','Gilad','ANTXR2','eSNP_DICE_ANTXR2_CD4_TH_data.csv'))
#divide eSNPs based on T-helper sub-population of effect
th1 = esnp_data.loc[esnp_data['Cell type'] == 'TH1']
th17 = esnp_data.loc[esnp_data['Cell type'] == 'TH17']
tfh = esnp_data.loc[esnp_data['Cell type'] == 'TFH']

import matplotlib.gridspec as gridspec

#Draw a barplot shwoint the effect-size of the ANTXR2 eSNPs
#on ANTXR2 expression per cell-type
fig = plt.figure(constrained_layout=True, figsize=(40, 100))
spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig)
b = fig.add_subplot(spec2[:1,:])
b.set_facecolor('w')
plt.bar(tfh['SNP ID'], tfh['Effect size'], width = 0.962, label='TFH', align = 'edge', linewidth = 0)
plt.bar(th17['SNP ID'], th17['Effect size'], width = 0.962, label='TH17', align = 'edge', linewidth = 0)
plt.bar(th1['SNP ID'], th1['Effect size'], width = 0.962, label='TH1', align = 'edge', linewidth = 0)
plt.xlim([0,len(tfh['SNP ID']) + len(th1['SNP ID'])])
plt.legend(loc='best', fontsize = 60)
plt.xticks(rotation='vertical', size = 20)
plt.yticks([-1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0],size = 60)
h = fig.add_subplot(spec2[2:,:])
ax = sns.heatmap(df_Rsq_to_heatmap, cmap = sns.color_palette("Reds"), square = True, linewidth = 0, cbar = False)
ax.tick_params(labelsize=60)
fig.savefig(os.path.join('C:\Users','barsh','Desktop','Gilad','ANTXR2','DICE_ANTXR2_Effect_Size_Barplot_and_ANTXR2_pairwise_LD_results_Rsq_heatmap.svg'))




