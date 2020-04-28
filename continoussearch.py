#!usr/bin/env python3
'''
Parse over all SSalignment files (individual runs and the overview).
Looking for the most common secondary  structure in the overview file.
Look for this secondary structure in all individual runs and seperate them (max_file).
The structure with the best energy (3D) is the one vor the next constrained run.

--force: instead of the most common secondary structure, the structure that fits
         the given constrain best

python3 continoussearch.py -p 00/surface/analyse/ --print --first 'CopStems_00.ss' --second 'CopStems_00_00_000000.ss' -f


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
    parser.add_argument ('--interaction',action='store_true', help='Find the interaction-structure most similar to the contrained one')
    parser.add_argument ('--first', help='First line in the dataframe') # CopStems_00.ss
    parser.add_argument ('--second', help='Second line in the dataframe') # CopStems_00_00_000000.ss
    parser.add_argument ('-i', '--initialname', help='CopStems_2rl_01, CopStems_2rl_02, ...')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level =  logging.DEBUG)
        print('DEBUGGING MODE')
    else:
        logging.basicConfig(level = logging.WARNING)

    path = args.path
    sumofruns = pd.DataFrame() #e.g CopStems_long-surface_10000.csv
    individual = pd.DataFrame() #e.g CopStems_long-surface_10000_01.csv

    individualpath = glob.glob(os.path.join(path, str(args.initialname+'_**.csv'))) #list of directory files
    individualpath.sort(key=lambda x: int(os.path.basename(x).split('.')[-2].split('_')[-1])) #nr of the file
    path_sumofruns = glob.glob(os.path.join(path, str(args.initialname+'.csv')))

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
    individual.drop(index_constrain, inplace=True)
    individual.drop(index_start, inplace=True)

    log.debug(sumofruns)
    log.debug(individual)

    if args.force:
        #get the most similar structure to the constraint
        count_constrain = sumofruns.min()['count_constraint'] #detect the min difference over all runs
        df_count_constrain = sumofruns[sumofruns['count_constraint']==sumofruns['count_constraint'].min()]

        print('count_constrain {}'.format(sumofruns.min()['count_constraint']))

        interim_best3d = {}

        for index, row in df_count_constrain.iterrows():
            interim_best3d[index] = row
            bp = interim_best3d[index][6]
            log.debug('bp {}'.format(bp))

            #Find in all runs the one with these bps and the minE
            df_constrain_sequence = (individual[individual['bp'] == bp]) #df with all these datasets from individual
            log.debug('df_constrain_sequence')
            log.debug(df_constrain_sequence)

        df_constrain_sequence.sort_values('energy_value',ascending=False) #sort energy
        log.debug('Full df_constrain_sequence')
        log.debug(df_constrain_sequence)

        best3d = (df_constrain_sequence.head(1))
        print('Best 3D structure: {}'.format(best3d))

    if args.interaction:
        #get the most similar structure to the interaction site
        count_interaction = sumofruns.min()['count_interaction_constraint'] #detect the min difference over all runs
        df_count_interaction = sumofruns[sumofruns['count_interaction_constraint']==sumofruns['count_interaction_constraint'].min()]

        #from this df geht the most similar structure to the whole constraint
        count_interaction_cc = df_count_interaction.min()['count_constraint'] #detect the min difference over all runs
        df_count_interaction_cc = df_count_interaction[df_count_interaction['count_constraint']==df_count_interaction['count_constraint'].min()]

        print('count_interaction {}, shape {} (rows, columns)'.format(count_interaction, df_count_interaction.shape))
        log.debug('df_count_interaction')
        log.debug(df_count_interaction)

        print('count_interaction_constraint {}, shape {} (rows, columns)'.format(count_interaction_cc, df_count_interaction_cc.shape))
        log.debug('df_count_interaction_cc')
        log.debug(df_count_interaction_cc)


        interim_best3d = {}
        # seperate for the 3D structure in the individual runs with the certain 2D structure
        for index, row in df_count_interaction_cc.iterrows():
            interim_best3d[index] = row
            bp = interim_best3d[index][6]
            #log.debug('bp {}'.format(bp))

            #Find in all runs the one with these bps and the minE
            df_interaction_cc_sequence = (individual[individual['bp'] == bp]) #df with all these datasets from individual
            #log.debug('df_interaction_cc_sequence')
            #log.debug(df_interaction_cc_sequence)

        df_interaction_cc_sequence.sort_values('energy_value',ascending=False) #sort energy
        log.debug('Full df_constrain_sequence')
        log.debug(df_interaction_cc_sequence)



        best3d = (df_interaction_cc_sequence.head(1))
        print('Best 3D structure: {}'.format(best3d))


    if not args.force and not args.interaction:
        #get the most common structure and search within the individual runs
        count_how_often = sumofruns.max()['count_how_often']
        sequence_of_count = sumofruns[sumofruns['count_how_often']==sumofruns['count_how_often'].max()]

        print('count {}'.format(sumofruns['count_how_often'].max()))

        interim_count = {}

        for index, row in sequence_of_count.iterrows():
            interim_count[index] = row
            bp = interim_count[index][6]
            log.debug('bp {}'.format(bp))

            df_max_sequence = (individual[individual['bp'] == bp])
            log.debug('df_max_sequence')
            log.debug(df_max_sequence)

        df_max_sequence.sort_values('energy_value',ascending=False)
        log.debug('Full df_max_sequence')
        log.debug(df_max_sequence)

        best3d = (df_max_sequence.head(1))
        print('Best 3D structure: {}'.format(best3d))

    #Extract the trafl-name and the list number from the .ss_detected name
    #CopStems_long-surface_10000_01-000109.ss_detected

    number = best3d.iloc[:,0].to_string(index=False).strip()
    log.debug('secondary structure file: {}'.format(number))
    filename = str((number.split('-')[0])+'-'+(number.split('-')[1])+'.trafl')
    line = int(number.split('.')[0].split('-')[-1])
    print('ss_file: {}, trafl_filename: {}, line: {}'.format(number,filename,line))
    bp_string = str(bp)

    file = open( str((number.split('-')[0])+'-'+(number.split('-')[1])+'.bp'), 'w')
    file.write(bp)
    file.close()

    if args.print:
        if args.force:
            df_constrain_sequence.to_csv(os.path.join(path,r'CopStems_constrain.csv'), mode='w',index=False, header=True,sep="\t")
        if args.interaction:
            df_interaction_cc_sequence.to_csv(os.path.join(path,r'CopStems_interaction_cc.csv'), mode='w',index=False, header=True,sep="\t")
        if not args.force and not args.interaction:
            df_max_sequence.to_csv(os.path.join(path,r'CopStems_common.csv'), mode='w',index=False, header=True,sep="\t")


if __name__ == "__main__":
    main()
