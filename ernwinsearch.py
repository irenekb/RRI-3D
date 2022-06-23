#!usr/bin/env python3

'''
Extract the structure from a ernwin log-file with the best minimum free
Energy.
Input:  out.log from an Ernwin run, with the following values:
        'Step','Sampling_Energy', 'Constituing_Energies','ROG','ACC','Asphericity',
        'Anisotropy','Local-Coverage','Tracked Energy','Tracked Energy',
        'Tracked Energy','time','Sampling Move','Rej.Clashes','Rej.BadMls'
Output: number of ernwin sample
'''

import logging
import argparse
import pandas as pd
#pd.set_option('display.max_columns', 500)
import csv
import glob
import os
import sys
import numpy as np
log = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description='Find minE file from a trafl-file')
    parser.add_argument ('-i', '--input', help='Path to Inputfiles')
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level =  logging.DEBUG)
        print('DEBUGGING MODE')
    else:
        logging.basicConfig(level = logging.WARNING)

    ernwin_df = pd.read_csv(args.input,sep='\t',index_col=0)

    #Ernwin seperator
    for col in ernwin_df.columns:
        if '|' in col:
            del ernwin_df[col]

    selernwin_df = ernwin_df.copy(deep=True)
    selernwin_df = selernwin_df[selernwin_df.index % 100==0]

    ernwin_df
    log.debug(ernwin_df)
    ernwin_df.sort_values('Sampling_Energy',ascending=True, inplace=True)
    bestenergy = (ernwin_df.head(1))
    step = ernwin_df[0:1].index.item()
    #step=bestenergy.iloc[0]['Step']
    log.debug('Best sampling energy:{}'.format(step))
    log.debug(bestenergy)

    selernwin_df
    log.debug(selernwin_df)
    selernwin_df.sort_values('Sampling_Energy',ascending=True, inplace=True)
    selbestenergy = (selernwin_df.head(1))
    selstep = selbestenergy[0:1].index.item()
    #selstep=selbestenergy.iloc[0]['Step']
    log.debug('Searching for sample: {}'.format(selstep))
    log.debug(selbestenergy)

    sys.exit(selstep)

if __name__ == "__main__":
    main()
