#!usr/bin/env python3

'''
Extract the structure from a traflfile with the best constrained minimum free
Energy.
Input .trafl an outputfile from SimRNA
Output min.trafl file
'''

import logging
import argparse
import csv
import pandas as pd
import glob
import os
log = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description='Find minE file from a trafl-file')
    parser.add_argument('-i', '--input', help='Input trafl-file')
    parser.add_argument ('-p', '--path', help='Path to Inputfiles')
    parser.add_argument ('-o', '--output', help='Output')
    parser.add_argument('-n', '--natural', action='store_true',
                        help='minE without constrained energy')
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level =  logging.DEBUG)
        print('DEBUGGING MODE')
    else:
        logging.basicConfig(level = logging.WARNING)

    path = args.path
    input = glob.glob(os.path.join(path, str(args.input)))
    output = glob.glob(os.path.join(path, str(args.output)))

    trafl_df = pd.DataFrame(columns=['consec_write_number','replica_number',
                                    'energy_value_plus_restraints_score','energy_value',
                                    'current_temperature','datapoints'])

    step = 2 #only need every second line
    with open(os.path.join(path,args.input)) as FILE:
        for nr, line in enumerate(FILE):
            if nr % step == 0:
                dat= [float(x) for x in line.split()]
                log.debug(dat)
            else:
                coord = line
                trafl_df = trafl_df.append(pd.Series([dat[0],dat[1],dat[2],
                                                    dat[3],dat[4],coord],
                                                    index=trafl_df.columns),
                                                    ignore_index=True)

    trafl_df.sort_values('energy_value',ascending=True, inplace=True)
    log.debug (trafl_df)

    bestenergy = (trafl_df.head(1))
    print('Best energy without constraint:')
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print (bestenergy)

    trafl_df.sort_values('energy_value_plus_restraints_score',ascending=True, inplace=True)
    bestenergy_cc = (trafl_df.head(1))
    print('Best energy with constraint:')
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(bestenergy_cc)

    coord = bestenergy_cc.iloc[0]['datapoints']
    bestenergy_cc = bestenergy_cc.drop('datapoints', 1)

    bestenergy_cc.to_csv(os.path.join(path,args.output), mode='w',index=False, header=False,sep=" ")

    with open(os.path.join(path,args.output), "a") as FILE:
        FILE.write(coord)

if __name__ == "__main__":
    main()
