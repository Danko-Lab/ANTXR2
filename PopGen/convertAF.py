"""
Convert VCF formatted .counts file to SF2 formatted AF file

"""

import csv
import os
import sys


def convert(file, path='/local/workdir/Iskander/ALLELE_COUNTS/'):

    #fileNAME = '{0}.chr4.phase3.frq.count'.format(file)
    fullPATH = os.path.join(path, '{0}.chr4.phase3.frq.count'.format(file))
    newFILE = '{0}.chr4.phase3.SF2.af'.format(file)
    with open(newFILE, 'w') as newAF:
        newAF.write('position\tx\tn\tfolded\n')
        with open(fullPATH, 'r') as vcfCounts:
            countsFILE = csv.reader(vcfCounts, delimiter='\t')
            next(countsFILE)

            for countsLINE in countsFILE:
                #Read in line from VCF counts file and then write it out
                polymorphism = countsLINE[5].split(':')[-1]#This will get the minor allele of the SNP
                newAF.write('{0}\t{1}\t{2}\t1\n'.format(countsLINE[1], polymorphism, countsLINE[3]))







        vcfCounts.close()
    newAF.close()

file=sys.argv[1]

convert(file)
