#!usr/bin/env python3
'''
Parse over all SSalignment files (individual runs and the overview).
Looking for the most common secondary  structure in the overview file.
Look for this secondary structure in all individual runs and seperate them (max_file).
The structure with the best energy (3D) is the one vor the next constrained run.

python3 continoussearch.py -p /scr/coridan/irene/Data_interaction/SimRNA_Start/CopStems_long-surface_10000/Analyse/ --print

To get the pdb for the next run:
SimRNA_trafl2pdbsstructure.pdb trajectory.trafl {list} AA
AA ... all atom reconstruction
list ... frame in the list
'''

import logging
import argparse
import csv
import os
import glob
import pandas as pd

log = logging.getLogger(__name__)
pd.options.display.max_colwidth = 500000000000 #otherwise it cut the strings


def main():
    parser = argparse.ArgumentParser(description='Find the best 3D structure after a SimRNA surface run and the SSAllignment analysis')
    parser.add_argument ('-p', '--path', help='Path to Inputfiles')
    parser.add_argument ('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument ('--print', action='store_true',help='Print a csv-file with all minEnergy relevant files')
    parser.add_argument ('-f','--force',action='store_true', help='Find the sec-structures most similar to the contrained one')
    parser.add_argument ('--first', help='First line in the dataframe') # CopStems_00.ss
    parser.add_argument ('--second', help='Second line in the dataframe') # CopStems_00_00_000000.ss
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level =  logging.DEBUG)
        print('DEBUGGING MODE')
    else:
        logging.basicConfig(level = logging.WARNING)

    path = args.path
    sumofruns = pd.DataFrame() #e.g CopStems_long-surface_10000.csv
    individual = pd.DataFrame() #e.g CopStems_long-surface_10000_01.csv

    individualpath = glob.glob(os.path.join(path, '*00_**.csv')) #list of directory files
    individualpath.sort(key=lambda x: int(os.path.basename(x).split('.')[-2].split('_')[-1])) #nr of the file
    path_sumofruns = glob.glob(os.path.join(path, '*_00.csv'))

    print(individualpath)
    print(path_sumofruns)

    #transfer to dataframe
    sumofruns = pd.read_csv(path_sumofruns[0],sep="\t",encoding='utf-8')

    for n, run in enumerate(individualpath):
        if n==0:
            individual = pd.read_csv(run,encoding='utf-8',sep="\t")

        if n!= 0:
            df = pd.read_csv(run,encoding='utf-8' ,sep="\t") #create dataframe for reading current csv
            individual = individual.append(df)#appends current csv to final DF

    # Get names of indexes for which column Age has value 30
    constrained_line = args.first
    start_line = args.second
    index_constrain = individual[individual['number'] == constrained_line].index
    index_start = individual[ individual['number'] == start_line].index

    #get the start and constrained ss and remove them from the df
    constrain_bp = individual.loc[individual['number'] == constrained_line].loc[0,'bp'].values[0]
    log.debug('constrain: {}'.format(constrain_bp))
    start_bp = individual.loc[individual['number'] == start_line].loc[1,'bp'].values[0]
    log.debug('start:  {}'.format(start_bp))
    # Delete these row indexes from dataFrame
    individual.drop(index_constrain , inplace=True)
    individual.drop(index_start , inplace=True)

    if args.force:
        df_count_constrain = sumofruns[sumofruns['count_constraint']==sumofruns['count_constraint'].min()]
        bp = df_count_constrain.iloc[:,6].to_string(index=False).strip() #bp with min dif to constrain
        #Find in all runs the one with these bps and the minE
        df_constrain_sequence = individual[individual['bp'] == bp] #df with all these datasets from individual
        df_constrain_sequence.sort_values('energy_value',ascending=False) #sort energy
        best_constrain = df_constrain_sequence.head(1) # take the one with the minE

        number = best_constrain.iloc[:,0].to_string(index=False).strip()
        filename = str((number.split('-')[0])+'-'+(number.split('-')[1])+'.trafl')
        #filename = str((number.split('-')[0])+'.trafl')
        line = int(number.split('.')[0].split('-')[2])

        print()
        print('count_constrain {}; Basepairs {}'.format(sumofruns.min()['count_constraint'],bp))
        print('ss_file: {}, trafl_filename: {}, line: {}'.format(number,filename,line))
        print()

        if args.print:
            df_constrain_sequence.to_csv(os.path.join(path,r'CopStems_long-surface_constraint.csv'), mode='w',index=False, header=True,sep="\t")

    log.debug(sumofruns)
    log.debug(individual)

    #get the most common structure and search within the individual runs
    count_how_often = sumofruns.max()['count_how_often']
    sequence_of_count = sumofruns[sumofruns['count_how_often']==sumofruns['count_how_often'].max()]
    bp = sequence_of_count.iloc[:,6].to_string(index=False).strip()

    log.debug(sequence_of_count)
    print('count {}; Basepairs {}'.format(count_how_often,bp))

    max_sequence = individual[individual['bp'] == bp]
    max_sequence.sort_values('energy_value',ascending=False) #sort energy
    best3d = max_sequence.head(1) # take the one with the minE

    log.debug(max_sequence)
    print(best3d)

    #Extract the trafl-name and the list number from the .ss_detected name
    #CopStems_long-surface_10000_01-000109.ss_detected

    number = best3d.iloc[:,0].to_string(index=False).strip()
    #filename = str((number.split('.')[0].split('-')[0])+'.csv')
    filename = str((number.split('-')[0])+'-'+(number.split('-')[1])+'.trafl')
    #filename = str((number.split('-')[0])+'.trafl')
    line = int(number.split('.')[0].split('-')[-1])
    print('ss_file: {}, trafl_filename: {}, line: {}'.format(number,filename,line))
    bp_string = str(bp)

    #df_bp = pd.DataFrame(bp)
    file = open( str((number.split('-')[0])+'-'+(number.split('-')[1])+'.bp'), 'w')
    file.write(bp)
    file.close()


    #get the most similar structure to the constraint

    #df1 = pandas.DataFrame(rawdat)
    if args.print:
        max_sequence.to_csv(os.path.join(path,r'CopStems_long-surface_max.csv'), mode='w',index=False, header=True,sep="\t")




if __name__ == "__main__":
    main()
