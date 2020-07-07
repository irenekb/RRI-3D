#!usr/bin/env python3

'''
    #python3 SimRNA-Scripts/trafls-analyse.py  -i /scr/coridan/irene/Data_interaction/SimRNA_Temp-constant/CopStems_inital_2buffer/output/CopStems_01_01.trafl -f CopStems_01 -t CopStems_01 -o n
'''
import logging
import argparse
import csv
from collections import OrderedDict
import os
import glob
import re
import copy
from operator import itemgetter
import itertools
from itertools import groupby, count
import sys
import pandas
from more_itertools import consecutive_groups
#from more-itertools import consecutive_groups

log = logging.getLogger(__name__)

def db2bps(db):
    """
    Get base pair list from db string.

    :param db: String of dotbracket
    """
    opening_bps = []
    bps = []

    for i in range(len(db)):
        if db[i] in ['(','[','{','<']:
            opening_bps.append(i)
        elif db[i] in [')',']','}','>']:
            bps.append([opening_bps.pop(),i]) # list of lists
    return bps

def trafl(traflfile):
    #fieldnames = ['row_in_trafl','consec_write_number', 'replica_number',
    #              'energy_values_plus_restraint_score', 'energy_value','current_temp']
    with open(traflfile, mode='r') as FILE:
        lines =  FILE.read().split('\n')
        coordinates = {}
        rawdat = {}

        line_count = 1
        for c, line in enumerate(lines):
            if not line:
                continue
            else:
                coor = []
                l = []
                entry = {}
                if (line_count % 2) == 0:
                    for value in line.split():
                        coor.append(float(value))
                    coordinates[c-1] = l
                    line_count += 1
                else:
                    for value in line.split():
                        l.append(float(value))
                    rawdat[int(l[0])] = l[2],l[3],l[4]
                    line_count += 1
    return rawdat

def difference(bigset,smalset,bp_dic):
    '''
    :param bigset: dictionary of eg all basepairs from the constraint, startpoint,
                  initial interaction, bp before...
    :param smalset: current set for comparison
    :param bp_dict: wich bps are unstable
    '''

    difference_bp= [x for x in smalset if x not in bigset] + [x for x in bigset if x not in smalset]
    count_difference = len(difference_bp)
    for pair in difference_bp:
        pair = str(pair)
        if pair in bp_dic.keys():
            bp_dic[pair] += 1
        else:
            bp_dic[pair] = 1

    return difference_bp, count_difference,bp_dic


def find_interaction(interim_interaction_db,chainbreak):
    '''
    Find the interacting basepairs

    :param interim_interaction_db: list of lists with the dotbracket structure
    :param chainbreak: int that marks the postion of the new chain
    :return interaction: list of lists with the basepairs
    '''

    findinteractionline = dict()

    for nr, db_line in enumerate(interim_interaction_db): #line includes only noncrossing bps
        bp_count = 0
        pairline = []
        interim_interaction_bps = db2bps(db_line)
        for pair in interim_interaction_bps: #counting interaction pairs
                if pair[0] < chainbreak and pair[1] > chainbreak:
                    bp_count += 1
        #al lines with an interaction with the bp count und the bps
        findinteractionline[nr] = [bp_count,interim_interaction_bps]
        log.debug(findinteractionline)

    #sort identified interactions ascending and take the last entry (most interaction bps)
    interaction_list = list(sorted(findinteractionline.values()))[-1]
    interaction = interaction_list[1]
    #how many bps are involved in the interaction
    interaction_countbp = interaction_list[0]
    #sort the interaction pairs ascending (first base)
    interaction.sort(key = lambda k: (k[0], -k[1]))

    #how long is the consecutive lenght of the interaction
    for pair in interaction:
        pairline.append(int(pair[0]))
        pairline.append(int(pair[1]))
    log.debug('pairline')
    log.debug(pairline)

    len_interaction = consecutive_interaction_length(pairline)

    #test_consecutive = [list(group)for group in consecutive_groups(pairline)]
    #log.debug('Consecutive Interaction')
    #log.debug(test_consecutive)
    log.debug('Interaction Length')
    log.debug(len_interaction)

    return interaction, interaction_countbp, len_interaction


def consecutive_interaction_length(pairline):
    '''
    Find the longest consecutive interaction

    :param pairline: list of all bases that are involved in the interaction
    :return len_interaction: length of the longest consectuive interaction
    '''
    pairset = list()
    c = count()

    for element in pairline: #Hash all array elements
        pairset.append(element)
    pairset.sort()
    log.debug('len pairline {}, pairline {}, pairset {}'.format(len(pairline),pairline,pairset))

    len_interaction = max((list(g) for _, g in groupby (pairset, lambda x : x-next(c))), key=len)
    log.debug('interaction len {}'.format(len_interaction))

    return (len_interaction)

def main():
    parser = argparse.ArgumentParser(description='Align SS from SimRNA')
    parser.add_argument ('-p', '--path', help='Path to Inputfiles')
    parser.add_argument ('-i', '--input', help='Input ss-sequence')
    parser.add_argument ('-c', '--constraint', help='Constrained ss-sequence')

    parser.add_argument ('-o', '--output', help='Name of the outputfile')
    parser.add_argument ('-m','--outputmode', choices=['w','a'], default='w', help="Overwrite ('w') or append ('a')")
    parser.add_argument ('-u', '--uniqueoutput', help='Name of the unique outputfile/or the already existing one')
    #parser.add_argument ('-b', '--bplist', help='Name of the bp-list (unique)/or the already existing one')

    parser.add_argument ('-t', '--trafl', help='Traflfile')

    parser.add_argument ('-v', '--verbose', action='store_true', help='Be verbose')

    #ToDo:
    parser.add_argument ('-r', '--right', action='store_true', help='expand right')
    parser.add_argument ('-l', '--left', action='store_true', help='expand_left')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level =  logging.DEBUG)
        print('DEBUGGING MODE')
    else:
        logging.basicConfig(level = logging.WARNING)

    output = args.output
    output_bp = output+'_bp'
    unique_outputfile = args.uniqueoutput
    unique_start_bp = args.uniqueoutput+'_bp'


    path = args.path
    start_ss = args.input
    constraint_ss = args.constraint
    start_fn = os.path.basename(args.input)


    dic = OrderedDict()
    dic_unique = OrderedDict() #dictionary for all already collected unique sequences with all the other data

    goal_start,goal_constraint,goal_interaction = 0,0,0 #How often do we get the contstrained secondary structure
    count_unique = {} #dictionary for bpstr and a counter
    unique_universal = {} #final dictonary for all unique sequences
    bp_dict_before = {} #collect the bps and how often do they change compared to the structure before
    bp_dict_start = {} #collect the bps and how often do they change compared to the start structure
    bp_dict_constraint = {} #collect the bps and how often do they change compared to the constrained structure
    bp_dict_interaction_cc = {} #collect the bps and how often do they change compared to the interaction itselfe


    #CONSTRAINED SECONDARY STRUCTURE
    with open(constraint_ss , 'r') as CSSFile:
        bp_constraint  = list()
        constraint_sequence = str()
        content = CSSFile.readlines()
        content = [s.rstrip() for s in content]

        interim_interaction_db =[]
        for line in content:
            bp_constraint.extend(db2bps(line))
            constraint_sequence = constraint_sequence +' '+ line
            interim_interaction_db.append(line)
            for nr, element in enumerate(line):
                if element == ' ':
                    log.debug('chainbreak: {}'.format(nr))
                    chainbreak = nr

        interaction_cc, interaction_countbp_cc, len_interaction_cc= find_interaction(interim_interaction_db,chainbreak)

        #bp_constraint = sorted(bp_constraint, key=lambda bp_constraint: bp_constraint[0])
        bp_constraint = sorted(bp_constraint, key=itemgetter(0))
        bp_constraint_str = ''.join(str(s) for s in bp_constraint)
        bp_interaction_str = ''.join(str(s) for s in interaction_cc)
        constraint_ssf = str(os.path.basename(constraint_ss))
        log.debug('constrain {}'.format(constraint_ssf))
        dic[constraint_ssf] = constraint_sequence,0,0,0,0,0,0,0,bp_constraint,0,0,0,0,interaction_cc,interaction_countbp_cc, len_interaction_cc,0,0
        # [0]         [1]               [2]         [3]         [4]       [5]           [6]       [7]      [8]
        # sequence,count_constraint,count_start,count_before,constancy,dif_constraint,dif_start,dif_before,bp

    #START SECONDARY STRUCTURE
    with open(start_ss , 'r') as SSFile:
        bp_start  = list()
        start_sequence = str()
        content = SSFile.readlines()
        content = [s.rstrip() for s in content]

        interim_interaction_db =[]
        for line in content:
            start_sequence = start_sequence  +' '+ line
            bp_start.extend(db2bps(line))
            interim_interaction_db.append(line)

        interaction, interaction_countbp, len_interaction = find_interaction(interim_interaction_db,chainbreak)

        #bp_start = sorted(bp_start, key=lambda bp_start: bp_start[0])
        bp_start = sorted(bp_start, key=itemgetter(0))
        bp_start_str = ''.join(str(s) for s in bp_start)
        bp_before = copy.deepcopy(bp_start)
        start_ssf = str(os.path.basename(start_ss))
        print('start {}'.format(start_ssf))

        dif_constraint,count_constraint,bp_dict_constraint = difference(bp_constraint,bp_start,bp_dict_constraint)

        dif_interaction_cc,count_interaction_constraint,bp_dict_interaction = difference(interaction_cc,interaction,bp_dict_interaction_cc)

        dic[start_ssf] = start_sequence, count_constraint,0,0,0,dif_constraint,0,0,bp_start,0,0,0,0,interaction, interaction_countbp, len_interaction, count_interaction_constraint,dif_interaction_cc
        # [0]         [1]               [2]         [3]         [4]       [5]           [6]       [7]      [8] [9]
        # sequence,count_constraint,count_start,count_before,constancy,dif_constraint,dif_start,dif_before,bp, time

    #APPENDMODE: if already values available - use mode append and take the values from their (unique)
    if args.outputmode == 'a':
        bp_dict_start = {}
        with open(args.uniqueoutput,'r') as Uniquefile:
            next(Uniquefile)
            for line in Uniquefile:
                seq,count_how_often,count_constraint,count_start,dif_constraint,dif_start,bp,bpstr,interaction,interaction_countbp,len_interaction,count_interaction_constraint,dif_interaction_cc = line.strip('\n').split('\t')
                count_unique[bpstr] = int(count_how_often)
                dic_unique[seq] = count_how_often,count_constraint,count_start, dif_constraint,dif_start,bp,bpstr,interaction,interaction_countbp,len_interaction,count_interaction_constraint,dif_interaction_cc
                #                       [0]                 [1]         [2]          [3]            [4]  [5] [6]   [7]              [8]              [9]                [10]                          [11]

        with open(unique_start_bp, 'r') as Bplist:
            next(Bplist)
            for line in Bplist:
                nr, bp1,bp2, count = line.strip('\n').split('\t')
                bp_dict_start[bp1,bp2] = int(count)

    print(path)
    print(start_ss)

    Files = glob.glob(os.path.join(path, '*.ss_detected')) #list of directory files
    Files.sort(key=lambda x: int(os.path.basename(x).split('.')[-2].split('-')[-1])) #nr of the file (for old and new)

    trafls = trafl(args.trafl)


    constancy = 0

    #START
    for File in Files:
        log.debug(File)
        with open(File,'r') as file:
            #get a correct bp List
            sequence = str()
            bp = list()
            content = file.readlines()
            content = [s.rstrip() for s in content]

            interim_interaction_db =[]
            for line in content:
                sequence = sequence +' '+ line
                bp.extend(db2bps(line))
                interim_interaction_db.append(line)
            bp = sorted(bp, key=itemgetter(0))

            filename = str(os.path.basename(File)) #exctract filename
            log.debug(filename)
            interaction,interaction_countbp, len_interaction = find_interaction(interim_interaction_db,chainbreak)

            #get the last entry from dict TODO: make it smarter
            last_key = str(list(dic.keys())[-1])
            last_entry = next(reversed(dic.values()))
            bp_before = last_entry[8] # basepairlist

            #If the current sequence differ to the one before safe it
            bpstr = ''.join(str(s) for s in bp)

            dif_constraint,count_constraint,bp_dict_interaction = difference(bp_constraint,bp,bp_dict_constraint)

            dif_interaction_cc, count_interaction_constraint,bp_dict_constraint = difference(interaction_cc,interaction,bp_dict_interaction_cc)

            dif_before, count_before, bp_dict_before = difference(bp_before,bp, bp_dict_before)

            dif_start, count_start, bp_dict_start = difference(bp_start,bp,bp_dict_start)


            if bp_start_str == bpstr: #check if we ever reach the constraint
                goal_start += 1

            if bp_constraint_str == bpstr: #check if we ever reach the constraint
                goal_constraint += 1

            if bp_interaction_str == bpstr: #check if we ever reach the interaction
                goal_interaction += 1

            #values for constancy - how long do we get the same sec.structure
            if count_before == 0:
                constancy += 1
            else:
                constancy =  0

            # get information from the corresponding trafl file
            timepoint_now = int(re.split('\W+',filename.split('-',-1)[1])[0])
            energy_values_plus_restraint_score = float()
            energy_value = float()
            current_temp = float()

            for key,value in trafls.items():
                if key == timepoint_now:
                    energy_values_plus_restraint_score = value[0]
                    energy_value = value[1]
                    current_temp = value[2]
                else:
                    continue

            # create the new entry in the dictonary
            dic[filename]=sequence, count_constraint, count_start, count_before,constancy,dif_constraint,dif_start,dif_before,bp,timepoint_now,energy_values_plus_restraint_score,energy_value,current_temp,interaction,interaction_countbp, len_interaction,count_interaction_constraint,dif_interaction_cc
            #                  [0]         [1]               [2]         [3]         [4]       [5]           [6]       [7]   [8]    [9]         [10]                                [11]           [12]     [13]            [14]                  [15]                       [16]

            log.debug('Difference to the constrained sequence, the start sequence {} and the one before {}'.format(count_constraint,count_start, count_before))

            #If the current sequence is an unique one safe it
            if bpstr not in count_unique.keys():
                count_unique[bpstr] = 1
                dic_unique[sequence] = 1,count_constraint,count_start,dif_constraint,dif_start,bp,bpstr,interaction,interaction_countbp,len_interaction,count_interaction_constraint,dif_interaction_cc
                #                       [0]                 [1]         [2]          [3]            [4]  [5] [6]   [7]              [8]              [9]                [10]                          [11]

            else: #increase the count for how often the structure appears
                for k,v in count_unique.items():
                    if k == bpstr:
                        count_unique[k] = v+1


    for k, v in dic_unique.items():
        for key, value in count_unique.items(): #[k]bpstr [v]count_how_often
            if key == v[6]:
                unique_universal[k] =value,v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11]

    #save the work
    with open(output, 'w') as csvfile:
        header = ['number', 'sequence','count_constraint',\
                  'count_start','count_before','constancy',\
                  'dif_constraint','dif_start','dif_before',\
                  'bp','time','energy_values_plus_restraint_score',\
                  'energy_value','current_temp','interaction',\
                  'interaction_countbp','len_interaction',\
                  'count_interaction_constraint','dif_interaction_cc']
        writer = csv.DictWriter(csvfile, fieldnames=header, delimiter = '\t')
        writer.writeheader()
        for k,v in dic.items():
            writer.writerow({'number':k,'sequence':v[0],'count_constraint':v[1],\
                             'count_start':v[2],'count_before':v[3],'constancy':v[4],\
                             'dif_constraint':v[5],'dif_start':v[6],'dif_before':v[7],\
                             'bp':v[8],'time':v[9],'energy_values_plus_restraint_score':v[10],\
                             'energy_value':v[11],'current_temp':v[12],'interaction':v[13],\
                             'interaction_countbp':v[14],'len_interaction': v[15],\
                             'count_interaction_constraint':v[16],'dif_interaction_cc':v[17]})

    #with open(unique_outputfile, args.outputmode) as csvfile:
    with open(unique_outputfile, 'w') as csvfile:
        header = ['sequence','count_how_often','count_constraint',\
                  'count_start','dif_constraint','dif_start',\
                  'bp','bpstr','interaction',\
                  'interaction_countbp','len_interaction',\
                  'count_interaction_constraint','dif_interaction_cc']
        writer = csv.DictWriter(csvfile, fieldnames=header, delimiter = '\t')
        writer.writeheader()
        for k,v in unique_universal.items():
            writer.writerow({'sequence':k,'count_how_often':v[0],'count_constraint':v[1],\
                             'count_start':v[2],'dif_constraint':v[3],'dif_start':v[4],\
                             'bp':v[5],'bpstr':v[6],'interaction':v[7],\
                             'interaction_countbp':v[8],'len_interaction':v[9],\
                             'count_interaction_constraint':v[10],'dif_interaction_cc':v[11]})

    header = ['nr','bp1','bp2', 'count']
    with open(output_bp, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header, delimiter = '\t')
        writer.writeheader()
        for k,v in enumerate(bp_dict_before.items()):
            writer.writerow({'nr':k,'bp1':v[0][1],'bp2':v[0][5],'count':v[1]})

    with open(unique_start_bp, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header, delimiter = '\t')
        writer.writeheader()
        for k,v in enumerate(bp_dict_before.items()):
            writer.writerow({'nr':k,'bp1':v[0][1],'bp2':v[0][5],'count':v[1]})

    print("We reach the constraint {} times.".format(goal_constraint))
    print("We reach the start {} times.".format(goal_start))
    print("We reach the interaction {} times.".format(goal_interaction))

if __name__ == "__main__":
    main()

'''NOTES:
Filename old
    CopStems_long-surface_10000_01-000001.ss_detected
    1. Which element
    2. surface from what
    3. Iterations
    4. Run
    5. Timepoint
    6. ss_detected

Filename new
    CopStems_expand_long00_03_1_1-005001.ss_detected
    1. Which element
    2. SimRNA script
    3. SimRNA config file
    4. interaction expansion round
    5. round for 'this' run
    6. seed
    7. Outputline
'''
