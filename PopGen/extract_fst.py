import os
import csv
import sys

def readP(inFile, path):
    """
    Read in .p files that are more extreme than 95% of distr and return their positions
    :param inFile:
    :param path:
    :return:
    """
    fullpath = os.path.join(path, inFile)
    SNP_pos = []
    with open(fullpath, 'r') as myP:
        csvP = csv.reader(myP, delimiter='\t')
        next(csvP)
        for field in csvP:
            if int(field[1]) <= 81060000 and int(field[1]) >= 81000000 and float(field[2]) <= .05:
                SNP_pos.append(int(field[1]))
                #print('chr4\t{0}\t{0}\t{1}'.format(int(field[1]), float(field[2])))
            else:
                pass
        myP.close()

    return SNP_pos

def regCalls(reg_file, path):
    '''
    This function will take in a BED file and retrieve all regulatory elements within our interesting region. There is no filtering for how confident the regulatory region is being called.
    That level of filtering can surely be easily implemented as another if/else statement or some such simple method.
    :param reg_file:
    :param path:
    :return:
    '''
    reg_intervals = []
    fullPath = os.path.join(path, reg_file)
    with open(fullPath, 'r') as myBED:
        for field in csv.reader(myBED, delimiter='\t'):
            if field[0] == 'chr4' and int(field[1]) >= 81000000 and int(field[2]) <= 81060000:#filter chr4 and our region of interest
                reg_intervals.append([int(field[1]), int(field[2])])
            else:
                pass
        myBED.close()

    return reg_intervals

def intersect_reg(reg_elements, fst_snp):
    called_snps = {}
    for snp in fst_snp:
        for element in reg_elements:
            if snp >= element[0] and snp <= element[1]:#if SNP is called in a regulatory region then we add it to our dictionary and break from loop
                called_snps[snp] = element# we assume that each SNP will hit only one regulatory region for the sake of efficiency and logically makes sense
                break
            else:
                pass
    return called_snps

def call_SNPs(snp, path, snp_file):
    snp_path = os.path.join(path, snp_file)
    with open(snp_path, 'r') as mySNP:
        snpRead = csv.reader(mySNP, delimiter='\t')
        next(snpRead)
        for field in snpRead:
            positions = int(field[1])
            try:#try statement to efficiently check if our SNP exists in hash if it does then add the allele frequencies
                snp[positions].append(field[4:])

            except KeyError:
                pass
        mySNP.close()
    return snp
def main():
    fst = ['CEU_YRI_snp.v2.p', 'CEU_CHB_snp.v2.p', 'CEU_JPT_snp.v2.p']#Modify this to change with fst.p files you want to extract high fst snps from
    local = '/home/iskander/Documents/Danko_lab/sweeps/FST_Estimates/P_VAL'
    fst_arrays = []
    for fst_file in fst:
        fst_positions = readP(inFile=fst_file, path=local)
        fst_arrays.append(set(fst_positions))
    fst_positions = list(set.intersection(fst_arrays[0], fst_arrays[1], fst_arrays[2]))

    reg_file = sys.argv[1]
    snp_file = sys.argv[2]


    reg_elements = regCalls(reg_file=reg_file, path='/local/storage/share-tmp/tc532_lac334/macs2_single/')
    called_snps = intersect_reg(reg_elements, fst_positions)
    snp_freqs = call_SNPs(snp=called_snps, path='/local/workdir/Iskander/ALLELE_COUNTS/', snp_file=snp_file)
    """
    Output formatting:
    chr4    snp_position    
    """
    candidate_snps = sorted(snp_freqs.keys())
    if len(candidate_snps) > 0:#handle no candidates found
        for snp in candidate_snps:
            snp_allele = ''
            for allele in snp_freqs[snp][2]:

                snp_allele = snp_allele + '\t{0}'.format(allele)

            print('chr4\t{0}{1}\t{2}\t{3}\t{4}'.format(snp, snp_allele, snp_freqs[snp][0], snp_freqs[snp][1], reg_file))
    else:
        pass




main()
