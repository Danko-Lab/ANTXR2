import os
import csv
import scipy.stats as stats
import sys


"""
Simple sript that will take a .windows.weir.fst formatted file from VCFtools and compute how extreme the mean Fst is in the distribution

"""

def plot_percentile(positions, fst):

    kde = stats.gaussian_kde(fst)
    for window in range(len(fst)):
        p = kde.integrate_box(fst[window], max(fst))
        print('{0}\t{1}'.format(positions[window], p))





def read(inFile, path):
    binStarts = []
    fst_arr = []
    fullPATH = os.path.join(path, inFile)
    with open(fullPATH, 'r') as myFst:
        readFst = csv.reader(myFst, delimiter='\t')
        next(readFst)
        for field in readFst:
            position = int(field[1])
            fst = float(field[5])
            binStarts.append(position)
            fst_arr.append(fst)

    myFst.close()

    return binStarts, fst_arr


def main():
    inF = sys.argv[1]
    inP = sys.argv[2]
    posits, fst = read(inF, inP)
    plot_percentile(posits, fst)

main()

