#!usr/bin/env python3

'''
Extract the structure from a traflfile with the best constrained minimum free
Energy.
Input .trafl an outputfile from SimRNA
Output min.trafl file

Values given by the trafl-file:
    'consec_write_number'
    'replica_number'
    'energy_value_plus_restraints_score'
    'energy_value'
    'current_temperature'
    'datapoints/cord' (newline)
 
'''

import logging
import argparse
import os
log = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description='Find minE file from a trafl-file')
    parser.add_argument('-i', '--input', help='Input trafl-file')
    parser.add_argument ('-p', '--path', help='Path to Inputfiles')
    parser.add_argument ('-o', '--output', help='Output')
    parser.add_argument('-n', '--natural', action='store_true',
                        help='minE without constrained energy') #TODO
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level =  logging.DEBUG)
        print('DEBUGGING MODE')
    else:
        logging.basicConfig(level = logging.WARNING)

    path = args.path
    cord_count = 1
    step = 2 #only need every second line

    with open(os.path.join(path,args.input)) as FILE:
        for nr, line in enumerate(FILE):
            if nr == 0:
                dat0= [float(x) for x in line.split()]

            elif nr % step == 0:
                dat1= [float(x) for x in line.split()]

                if dat1[3] < dat0[3]:
                    dat0 = [float(x) for x in line.split()]
                    cord_count = nr+1

            elif nr == cord_count:
                cord = line

            else:
                continue
    
    with open (os.path.join(path,args.output), "w") as FILE:
        for element in dat0:
            FILE.write("{} ".format(element))
        FILE.write("\n")
        FILE.write("{}".format(cord))

if __name__ == "__main__":
    main()
