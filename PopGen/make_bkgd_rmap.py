import csv

import argparse

"""
------------------------------------------------------------------------------------
MIT License

Copyright (c) 2018 Iskander Said

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------------

A modular command line parser for the coder on the go. 

Simply add or remove arguments depending on your program. The code implementation is quite straightforward simply
write:

def main():
#Call command line class
    myCommandLine = CommandLine()
    #Call a commandline argument
    myArg = myCommandLine.args.myArg

main()

Done. Boom. 
"""


class CommandLine():
    def __init__(self, inOpts=None):
        self.parser = argparse.ArgumentParser(
            description='Simple script to rewrite bed files into rmap or bkgd files just make sure the info field in the bed files is the recomb rate, or background selection value',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            # Allows the epilog to be formatted in the way I want!!
            epilog=('''                                      
              Rewrite rmap or bkgd bed files

              '''),
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
        )  ###Add or remove optional arguments here
        self.parser.add_argument("-r", "--rmap", type=str, action="store", nargs="?", help="Recombination map file",
                                 default=False)
        self.parser.add_argument("-bkgd", "--bkgd", type=str, action="store", nargs="?", help="Background selection file",
                                 default=False)
        self.parser.add_argument("-af", "--allele_freq", type=str, action="store", nargs="?",
                                 help="Allele freq file",
                                 default=False)
        self.parser.add_argument("-o", "--out", type=str, action="store", nargs="?",
                                 help="Output file name",
                                 default='output')
        if inOpts is None:  ## DONT REMOVE THIS. I DONT KNOW WHAT IT DOES BUT IT BREAKS IT
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


def parseB(inFile, bkgd, output):
    #convert bed back bkgd
    bHASH = {}  # Read ranges into memory
    with open(bkgd, 'r') as myB:
        for field in csv.reader(myB, delimiter='\t'):
            bHASH[(int(field[1]), int(field[2]))] = float(field[3])
        myB.close()

    with open('{0}.bkgd'.format(output), 'w') as myBkg:
        myBkg.write('position\tbvalue\n')
        with open('{0}'.format(inFile), 'r') as myAF:  # Iterate through values
            AF_file = csv.reader(myAF, delimiter='\t')
            next(AF_file)
            for line in AF_file:
                position = int(line[0])
                for ranges in bHASH.keys():
                    if position < ranges[1] and position > ranges[0]:
                        myBkg.write('{0}\t{1}\n'.format(position, bHASH[(ranges[0], ranges[1])]))
                        break
                    else:
                        pass
            myAF.close()
        myBkg.close()



def parseRMAP(inFile, rmap, output):
    #convert bed back to rmap
    bHASH = {}  # Read ranges into memory
    with open(rmap, 'r') as myB:
        for field in csv.reader(myB, delimiter='\t'):
            bHASH[(int(field[1]), int(field[2]))] = float(field[3])
        myB.close()

    with open('{0}.rmap'.format(output), 'w') as myBkg:
        myBkg.write('position\trate\n')
        with open('{0}'.format(inFile), 'r') as myAF:  # Iterate through values
            AF_file = csv.reader(myAF, delimiter='\t')
            next(AF_file)
            for line in AF_file:
                position = int(line[0])
                for ranges in bHASH.keys():
                    if position < ranges[1] and position > ranges[0]:
                        myBkg.write('{0}\t{1}\n'.format(position, bHASH[(ranges[0], ranges[1])]))
                        break
                    else:
                        pass
            myAF.close()
        myBkg.close()



def main():
    myCommandLine = CommandLine()
    if myCommandLine.args.rmap != False:
        parseRMAP(inFile=myCommandLine.args.allele_freq, rmap=myCommandLine.args.rmap, output=myCommandLine.args.out)
    elif myCommandLine.args.bkgd != False:
        parseB(inFile=myCommandLine.args.allele_freq, bkgd=myCommandLine.args.bkgd, output=myCommandLine.args.out)
main()