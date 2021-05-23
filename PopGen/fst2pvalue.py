# This script will compute the Gaussian KDE of the Fst on the chromosome and then determine the p-value of every SNP in the file
# Then we will output the p-values in file readable by QQPLOT
import os
import csv
import scipy.stats as stats
import sys


def read_file(inFile, path):  # Read in the .weir.fst file
    fullPATH = os.path.join(path, inFile)
    snp_position = []
    fstArray = []
    with open(fullPATH, 'r') as myFile:
        fst = csv.reader(myFile, delimiter='\t')
        next(fst)
        for field in fst:
            fstArray.append(max(float(field[2]), 0))
            snp_position.append(int(field[1]))
        myFile.close()
    fst_snp = [fstArray, snp_position]

    return fst_snp


def kde_function(fst_values):  # Create a gaussian kernel density estimate PDF from which to calculate the likelihoods

    kde = stats.gaussian_kde(fst_values[0])

    return kde


def compute_oneSided_pvalue(pdf, point, tail):
    return pdf.integrate_box(point, tail)


def writeFile(pvalues, positions, file_desig):
    with open(file_desig, 'w') as myP:
        myP.write('CHROM\tPOS\tPVALUE\n')
        for p in range(len(pvalues)):
            myP.write('4\t{0}\t{1}\n'.format(positions[p], pvalues[p]))
        myP.close()


def main(inFile, path):
    pval_array = []
    fst_snp = read_file(inFile, path)
    gauss = kde_function(fst_snp)
    tail_end = max(fst_snp[0])
    for fst in fst_snp[0]:
        pvalue = compute_oneSided_pvalue(pdf=gauss, point=fst, tail=tail_end)
        pval_array.append(pvalue)
    new_file = inFile.split('.')[0] + '.p'
    writeFile(pvalues=pval_array, positions=fst_snp[1], file_desig=new_file)

input_file = sys.argv[1]
input_path = sys.argv[2]

main(input_file, input_path)