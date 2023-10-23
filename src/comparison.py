#!usr/bin/env python3

import logging
import argparse
import csv
from collections import OrderedDict
import os
import glob
import re
import copy
from operator import itemgetter
#import itertools
import pandas as pd
from more_itertools import consecutive_groups

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
    """
    Parse the SimRNA trajectory file

    :param traflfile: *.trafl file with the fields: ['row_in_trafl',
                        'consec_write_number', 'replica_number',
                        'energy_values_plus_restraint_score',
                        'energy_value','current_temp']
    """

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

def find_difference(bigset,smalset):
    '''
    :param bigset: dictionary of eg all basepairs from the constraint, startpoint,
                  initial interaction, bp before...
    :param smalset: current set for comparison
    '''
    log.debug("bigset {}".format(bigset))
    log.debug("smalset {}".format(smalset))

    if smalset == 0:
        difference_bp = bigset
        count_difference = len(bigset)
    else:
        difference_bp= [x for x in smalset if x not in bigset] + [x for x in bigset if x not in smalset]
        count_difference = len(difference_bp)

    return difference_bp, count_difference


def find_interaction(interim_interaction_db,chainbreak):
    '''
    Find the interacting basepairs

    :param interim_interaction_db: list of lists with the dotbracket structure
    :param chainbreak: int that marks the postion of the new chain
    :return interaction: list of lists with the basepairs
    '''

    intra = list()
    consec_sum = list()

    nonconsecutive = dict()

    for nr, db_line in enumerate(interim_interaction_db): #line includes only noncrossing bps
        log.debug(nr)
        log.debug(db_line)

        c =list()
        consecutive = list()

        interim_interaction_bps = db2bps(db_line)
        interim_interaction_bps.sort()

        stop = len(interim_interaction_bps)-1

        count1 = 0
        for count2, pair in enumerate(interim_interaction_bps):
            log.debug('count1 {} count2 {} pair {} cb {}'.format(count1, count2, pair,chainbreak))

            #counting interaction pairs
            if pair[0] < chainbreak and pair[1] > chainbreak:
                if count1 == 0:
                    bp0,bp1 = pair[0],pair[1]
                    c.append(pair)
                    count1 += 1
                    #print('start {} {}'.format(bp0,bp1))
                else:
                    if pair[0] == bp0+1 and pair[1] == bp1-1:
                        #print('append {} {}'.format(pair[0],pair[1]))
                        c.append(pair)
                        if count2 == stop:
                            consecutive.append(c.copy())
                            consec_sum.append(c.copy())
                            c.clear()

                    else:
                        consecutive.append(c.copy())
                        consec_sum.append(c.copy())
                        c.clear()
                        c.append(pair)

                    bp0,bp1 = pair[0],pair[1]


            #collect all intramolecular bases
            else:
                intra.append(pair[0])
                intra.append(pair[1])

            if count2 == stop:
                if c:
                    consecutive.append(c.copy())
                    consec_sum.append(c.copy())
                    c.clear()

        log.debug('consecutive stems {}'.format(consecutive))

        if consecutive:
            pos = 1
            for conseccount, cstem in enumerate(consecutive):
                log.debug('stem  to check {} {}'.format(conseccount,cstem))
                if conseccount == 0:
                    pos = 1
                    stem1 = cstem
                    start1 = stem1[-1][0]
                    end1 = stem1[-1][1]
                    nonconsecutive[(str(nr) + str(pos))] = [stem1]
                    log.debug('nr {} , start1 {}, end1 {}'.format((str(nr) + str(pos)), start1, end1))

                else:
                    stem2 = cstem
                    start2 = stem2[0][0]
                    end2 = stem2[0][1]
                    log.debug('nr {} , start1 {}, end1 {}, start2 {}, end2 {}'.format((str(nr) + str(pos)), start1, end1, start2, end2))
                    check =list()

                    for i in intra:
                        if (i > start1 and i < start2) or (i > end2 and i < end1):
                            check.append(i)
                    if not check:
                        nonconsecutive[(str(nr) + str(pos))].append(stem2)
                    else:

                        pos += 1
                        nonconsecutive[(str(nr) + str(pos))] = [stem2]
                    start1 = stem2[-1][0]
                    end1 = stem2[-1][1]

    log.debug('all consec {} '.format(consec_sum) )
    for k, ss in nonconsecutive.items():
        log.debug('nonconsecutive pos {} {}'.format(k, ss))
    log.debug('intras {}'.format(intra))

    if consec_sum:
        maxconsecutive_bp = max((x) for x in consec_sum)
        maxconsecutive_len = max(len(x) for x in consec_sum)
    else:
        maxconsecutive_bp = 0
        maxconsecutive_len = 0
    length = 0

    if nonconsecutive:
        for k, currentbps in nonconsecutive.items():
            currentbps = [j for i in currentbps for j in i]
            if len(currentbps) > length:
                maxinteraction_bp = currentbps
                maxinteraction_len = len(currentbps)
                length = maxinteraction_len
    else:
        maxinteraction_bp = 0
        maxinteraction_len = 0

    log.debug('RESULT')
    log.debug(maxinteraction_bp,  maxinteraction_len, maxconsecutive_len,maxconsecutive_bp)
    return maxinteraction_bp,  maxinteraction_len, maxconsecutive_len
    #return interaction, interaction_countbp, len_interaction


def find_intra(bp,interaction_countbp,chainbreak):
    '''
    Return the intramolecular basepairs for both chains seperatly
    '''

    bp_set = set(map(tuple, bp))

    if interaction_countbp == 0:
        bp_rest = bp_set
    else:
        interaction_countbp_set = set(map(tuple, interaction_countbp))
        #remove already identified interaction pairs
        bp_rest = bp_set.symmetric_difference(interaction_countbp_set)

    intra_chainA = list()
    intra_chainB = list()
    otherinter = list()

    for bp in bp_rest:
        if bp[0] < chainbreak and bp[1] < chainbreak:
            intra_chainA.append(bp)
        if bp[0] > chainbreak and bp[1] > chainbreak:
            intra_chainB.append(bp)
        if bp[0] < chainbreak and bp[1] > chainbreak:
            otherinter.append(bp)

    return intra_chainA, len(intra_chainA), intra_chainB, len(intra_chainB), otherinter, len(otherinter)


def find_not_in_path(bp,set_constrained_se_bps):
    '''
    Return basepairs that are not listed in the start or the end/target    

    param bp: bp list of the current structure
    param set_constrained_se_bps: bp set with all bps from the initial start and the end/target structure
    '''

    bp = set(map(tuple, bp))
    not_in_path = list(bp.difference(set_constrained_se_bps))

    return not_in_path, len(not_in_path)

def find_multiplets(bps):
    '''
    Return a list of multiplets

    param bp: bp list of the current structure
    '''
    flat_list = [item for sublist in bps for item in sublist]
    duplicates = list(set([x for x in flat_list if flat_list.count(x) > 1]))

    multiplet = list()

    for bp in bps:
        for d in duplicates:
            if bp[0] == d or bp[1] == d:
                if bp not in  multiplet:
                    multiplet.append(bp)
            else:
                continue

    return multiplet


def main():
    parser = argparse.ArgumentParser(description='Align SS from SimRNA')
    parser.add_argument ('-p', '--path', help='Path to Inputfiles')
    parser.add_argument ('-i', '--input', help='Input ss-sequence')
    parser.add_argument ('-c', '--constraint', help='Constrained ss-sequence')
    parser.add_argument ('-o', '--output', help='Name of the outputfile')
    parser.add_argument ('-m','--outputmode', choices=['w','a'], default='w',
                        help="Overwrite ('w') or append ('a')")
    parser.add_argument ('-u', '--uniqueoutput',
                        help='Name of the unique outputfile/or the already existing one')
    parser.add_argument ('-t', '--trafl', help='Traflfile')
    parser.add_argument ('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument ('-s', '--start', help='Path to the inital starting ss-file')
    parser.add_argument ('-e', '--end', help='Path to the end/target ss-file')

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

    dic = OrderedDict()
    dic_unique = OrderedDict() #dictionary for all already collected unique sequences with all the other data

    #COUNTER: How often do we get the contstrained secondary structure
    goal_start,goal_constraint,goal_interaction = 0,0,0

    count_unique = {} #dictionary for bpstr and a counter
    unique_universal = {} #final dictonary for all unique sequences

    bp_dict_currentrun = {}
    bp_dict_unique = {}

    #INITIAL START STRUCTURE
    with open(args.start , 'r') as STARTFILE:
        bp_initstart  = list()
        content = STARTFILE.readlines()
        content = [s.rstrip() for s in content]

        for line in content:
            bp_initstart.extend(db2bps(line))
    #start_set = set(map(tuple, bp_start))

    #END/TARGET STRUCTURE
    with open(args.end , 'r') as ENDFILE:
        bp_end  = list()
        content = ENDFILE.readlines()
        content = [s.rstrip() for s in content]

        for line in content:
            bp_end.extend(db2bps(line))


    print(bp_initstart)
    print(bp_end)
    #prep. for not in path function
    together = bp_initstart + bp_end
    constrained_se_bps = [] 
    [constrained_se_bps.append(x) for x in together if x not in constrained_se_bps] 
    set_constrained_se_bps = set(map(tuple, constrained_se_bps))

    print(set_constrained_se_bps)

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
            #severalchains=False #only for non-interaction calculations
            for nr, element in enumerate(line):
                if element == ' ':
                    chainbreak = nr
                    severalchains=True
                    log.debug('chainbreak: {}'.format(nr))
            #TODO:debug
            #if severalchains == False:
                #chainbreak = len(line)+1
                #print('Only one chain available!')

        interaction_cc, interaction_countbp_cc, len_interaction_cc= find_interaction(interim_interaction_db,chainbreak)

        #if args.intra: 
        intra_chainA,  intra_chainA_count, intra_chainB, intra_chainB_count, other_inter, other_inter_count= find_intra(bp_constraint,interaction_cc,chainbreak)
        not_in_path, not_in_path_count = find_not_in_path(bp_constraint, set_constrained_se_bps)
        multiplet = find_multiplets(bp_constraint)

        #bp_constraint = sorted(bp_constraint, key=lambda bp_constraint: bp_constraint[0])
        bp_constraint = sorted(bp_constraint, key=itemgetter(0))
        bp_constraint_str = ''.join(str(s) for s in bp_constraint)
        bp_interaction_str = ''.join(str(s) for s in interaction_cc)
        constraint_ssf = str(os.path.basename(constraint_ss))
        log.debug('constrain {}'.format(constraint_ssf))
        dic[constraint_ssf] = constraint_sequence,0,0,0,0,0,0,0,bp_constraint,0,0,0,0,interaction_cc,interaction_countbp_cc, len_interaction_cc,0,0,intra_chainA,  intra_chainA_count, intra_chainB, intra_chainB_count, other_inter, other_inter_count,not_in_path, not_in_path_count,multiplet
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
        intra_chainA,  intra_chainA_count, intra_chainB, intra_chainB_count, other_inter, other_inter_count = find_intra(bp_start,interaction,chainbreak)
        not_in_path, not_in_path_count = find_not_in_path(bp_start, set_constrained_se_bps)
        multiplet = find_multiplets(bp_start)

        #bp_start = sorted(bp_start, key=lambda bp_start: bp_start[0])
        bp_start = sorted(bp_start, key=itemgetter(0))
        bp_start_str = ''.join(str(s) for s in bp_start)
        bp_before = copy.deepcopy(bp_start)
        start_ssf = str(os.path.basename(start_ss))
        print('start {}'.format(start_ssf))

        dif_constraint,count_constraint = find_difference(bp_constraint,bp_start)

        dif_interaction_cc,count_interaction_constraint = find_difference(interaction_cc,interaction)

        dic[start_ssf] = start_sequence, count_constraint,0,0,0,dif_constraint,0,0,bp_start,0,0,0,0,interaction, interaction_countbp, len_interaction, count_interaction_constraint,dif_interaction_cc,intra_chainA,  intra_chainA_count, intra_chainB, intra_chainB_count, other_inter, other_inter_count,not_in_path, not_in_path_count,multiplet
        # [0]         [1]               [2]         [3]         [4]       [5]           [6]       [7]      [8] [9]
        # sequence,count_constraint,count_start,count_before,constancy,dif_constraint,dif_start,dif_before,bp, time

    #APPENDMODE: if already values available - use mode append and take the values from their (unique)
    if args.outputmode == 'a':
        bp_dict_unique = {}
        with open(args.uniqueoutput,'r') as Uniquefile:
            next(Uniquefile)
            for line in Uniquefile:
                seq,count_how_often,count_constraint,count_start,dif_constraint,dif_start,bp,bpstr,interaction,interaction_countbp,len_interaction,count_interaction_constraint,dif_interaction_cc,intra_chainA,  intra_chainA_count, intra_chainB, intra_chainB_count, other_inter, other_inter_count,not_in_path, not_in_path_count,multiplet = line.strip('\n').split('\t')
                count_unique[bpstr] = int(count_how_often)
                dic_unique[seq] = count_how_often,count_constraint,count_start, dif_constraint,dif_start,bp,bpstr,interaction,interaction_countbp,len_interaction,count_interaction_constraint,dif_interaction_cc,intra_chainA,  intra_chainA_count, intra_chainB, intra_chainB_count, other_inter, other_inter_count, not_in_path, not_in_path_count,multiplet
                #                       [0]                 [1]         [2]          [3]            [4]  [5] [6]   [7]              [8]              [9]                [10]                          [11]

        with open(unique_start_bp, 'r') as Bplist:
            next(Bplist)
            for line in Bplist:
                nr, bp1,bp2, count = line.strip('\n').split('\t')
                pairstr = str(bp1+', '+bp2)
                bp_dict_unique[pairstr] = int(count)


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
            intra_chainA,  intra_chainA_count, intra_chainB, intra_chainB_count, other_inter, other_inter_count= find_intra(bp,interaction,chainbreak)
            not_in_path, not_in_path_count = find_not_in_path(bp, set_constrained_se_bps)
            multiplet = find_multiplets(bp)

            #get the last entry from dict TODO: make it smarter
            last_key = str(list(dic.keys())[-1])
            last_entry = next(reversed(dic.values()))
            bp_before = last_entry[8] # basepairlist

            #If the current sequence differ to the one before safe it
            bpstr = ''.join(str(s) for s in bp)

            log.debug('check 1 {} {}'.format(bp_constraint,bp))

            dif_constraint,count_constraint = find_difference(bp_constraint,bp)

            log.debug('check 2 {} {}'.format(interaction_cc,interaction))

            dif_interaction_cc, count_interaction_constraint = find_difference(interaction_cc,interaction)
            dif_before, count_before = find_difference(bp_before,bp)
            dif_start, count_start = find_difference(bp_start,bp)


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
            dic[filename]=sequence, count_constraint, count_start, count_before,constancy,dif_constraint,dif_start,dif_before,bp,timepoint_now,energy_values_plus_restraint_score,energy_value,current_temp,interaction,interaction_countbp, len_interaction,count_interaction_constraint,dif_interaction_cc,intra_chainA, intra_chainA_count, intra_chainB, intra_chainB_count, other_inter, other_inter_count,not_in_path, not_in_path_count,multiplet
            #                  [0]         [1]               [2]         [3]         [4]       [5]           [6]       [7]   [8]    [9]         [10]                                [11]           [12]     [13]            [14]                  [15]                       [16]

            bp = sorted(bp, key=itemgetter(0))
            for pair in bp:
                pairstr=(', '.join(str(x) for x in pair))
                bp_dict_currentrun[pairstr] = bp_dict_currentrun.get(pairstr, 0) + 1

            for k, v in bp_dict_currentrun.items():
                log.debug('bp_dict_currentrun: k {}, v {}'.format(k,v))
                if k in bp_dict_unique:
                    log.debug('adding: k {} : +v {}'.format(k,v))
                    bp_dict_unique[k] = v+1
                else:
                    log.debug('new: k {} : +v {}'.format(k,v))
                    bp_dict_unique[k] = v


            log.debug('Difference to the constrained sequence, the start sequence {} and the one before {}'.format(count_constraint,count_start, count_before))

            #If the current sequence is an unique one safe it
            if bpstr not in count_unique.keys():
                count_unique[bpstr] = 1
                dic_unique[sequence] = 1,count_constraint,count_start,dif_constraint,dif_start,bp,bpstr,interaction,interaction_countbp,len_interaction,count_interaction_constraint,dif_interaction_cc,intra_chainA,  intra_chainA_count, intra_chainB, intra_chainB_count, other_inter, other_inter_count,not_in_path, not_in_path_count,multiplet
                #                       [0]                 [1]         [2]          [3]            [4]  [5] [6]   [7]              [8]              [9]                [10]                          [11]
            else: #increase the count for how often the structure appears
                for k,v in count_unique.items():
                    if k == bpstr:
                        count_unique[k] = v+1


    for k, v in dic_unique.items():
        for key, value in count_unique.items(): #[k]bpstr [v]count_how_often
            if key == v[6]:
                unique_universal[k] =value,v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15],v[16],v[17],v[18],v[19],v[20]


    columnnames_individual = ['number','sequence','count_constraint',
                              'count_start','count_before','constancy',
                              'dif_constraint','dif_start','dif_before',
                              'bp','time','energy_values_plus_restraint_score',
                              'energy_value','current_temp','interaction',
                              'interaction_countbp','len_interaction',
                              'count_interaction_constraint','dif_interaction_cc',
                              'intra_chainA',  'intra_chainA_count', 
                              'intra_chainB', 'intra_chainB_count',
                              'other_inter', 'other_inter_count',
                              'not_in_path', 'not_in_path_count',
                              'multiplet']

    columnnames_unique = ['sequence','count_how_often','count_constraint',
                          'count_start','dif_constraint','dif_start',
                          'bp','bpstr','interaction',
                          'interaction_countbp','len_interaction',
                          'count_interaction_constraint','dif_interaction_cc',
                          'intra_chainA',  'intra_chainA_count', 
                          'intra_chainB', 'intra_chainB_count',
                          'other_inter', 'other_inter_count',
                          'not_in_path', 'not_in_path_count',
                          'multiplet']


    df_individual = pd.DataFrame.from_dict(dic,orient='index').reset_index()
    df_individual.columns = columnnames_individual
    df_unique = pd.DataFrame.from_dict(unique_universal, orient='index' ).reset_index()
    df_unique.columns = columnnames_unique


    # save the files
    df_individual.to_csv(output, sep='\t',index=False)

    df_unique.to_csv(unique_outputfile, sep='\t',index=False)

    header = ['nr','bp1','bp2', 'count']
    with open(output_bp, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header, delimiter = '\t')
        writer.writeheader()
        for k,v in enumerate(bp_dict_currentrun.items()):
            currentpair=[x.strip() for x in v[0].split(', ')]
            writer.writerow({'nr':k,'bp1':currentpair[0],'bp2':currentpair[1],'count':v[1]})

    with open(unique_start_bp, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header, delimiter = '\t')
        writer.writeheader()
        for k,v in enumerate(bp_dict_unique.items()):
            currentpair= currentpair=[x.strip() for x in v[0].split(', ')]
            writer.writerow({'nr':k,'bp1':currentpair[0],'bp2':currentpair[1],'count':v[1]})

    print("We reach the constraint {} times.".format(goal_constraint))
    print("We reach the start {} times.".format(goal_start))
    print("We reach the interaction {} times.".format(goal_interaction))

if __name__ == "__main__":
    main()
