import csv
import sys

def intersectFiles(bkgd, rmap, af, output):
    newAF = af.split('.')[0] + '.pruned.af'
    bHASH = {}
    with open(bkgd, 'r') as myB:  # Create a set to find intersections and create a hash to store the bvalues
        for field in csv.reader(myB, delimiter='\t'):
            if field[0] not in bHASH.keys():
                bHASH[int(field[0])] = field[1]
        myB.close()

    afDict = {}
    with open(af, 'r') as myAF:  # To remove the extra .af fields
        # We will use keys to create sets and find the intersections here; this hash will be able to track duplicates
        AFreader = csv.reader(myAF, delimiter='\t')
        next(AFreader)
        for line in AFreader:
            if int(line[0]) not in afDict.keys():
                afDict[int(line[0])] = [[line[1], line[2], line[3]]]  # in case where there are no duplicates
            else:
                afDict[int(line[0])].append([line[1], line[2], line[3]])  # in case where there are duplicate positions
        myAF.close()

    rHASH = {}
    with open(rmap, 'r') as myRMAP:
        for field in csv.reader(myRMAP, delimiter='\t'):
            rHASH[int(field[0])] = field[1]
        myRMAP.close()

    # We have loaded all of these hash tables and now lets create the intersection
    intersection = set(afDict.keys()) & set(bHASH.keys()) & set(rHASH.keys())
    intersection = sorted(list(intersection))
    with open('{0}.rmap'.format(output), 'w') as myRm:
        myRm.write('position\trate\n')
        with open('{0}.bkgd'.format(output), 'w') as myBk:
            myBk.write('position\tbvalue\n')
            with open(newAF, 'w') as myAF:
                myAF.write('position\tx\tn\tfolded\n')
                for position in intersection:
                    for field in afDict[position]:
                        myAF.write('{0}\t{1}\t{2}\t{3}\n'.format(position, field[0], field[1], field[2]))
                        myBk.write('{0}\t{1}\n'.format(position, bHASH[position]))
                        myRm.write('{0}\t{1}\n'.format(position, rHASH[position]))
                myAF.close()
            myBk.close()
        myRm.close()



def main():
    bkgd = sys.argv[1]
    rmap = sys.argv[2]
    af = sys.argv[3]
    out = sys.argv[4]

    intersectFiles(bkgd, rmap, af, out)
