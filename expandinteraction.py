#!usr/bin/env python3
'''
TESTINPUT
python3 expandinteraction.py  -b 2 -n test5.seq -x test5.bp -o test5_01.ss
python3 expandinteraction.py  -b 2 -n test5.seq -d test5.ss -o test5_02.ss

python3 expandinteraction.py  -b 2 -n SimRNA_interaction/1zci.seq -x SimRNA_interaction/1zci.bp -o 1zci_01.ss
python3 expandinteraction.py  -b 2 -n SimRNA_interaction/1zci.seq -d SimRNA_interaction/1zci.ss -o 1zci_02.ss
'''

import logging
import argparse
import csv
from collections import defaultdict
import json
from operator import itemgetter
import sys
import numpy as np
from sys import exit
np.set_printoptions(threshold=sys.maxsize,linewidth =900)


log = logging.getLogger(__name__)

def is_crossing(dotbracket,basepair):
    """
    Check if the current baspair is crossing or not (True|False)

    :param dotbracket: Current dotbracketline (list)
    :param baspair: Current basepair
    """
    count = 0
    crossing = True

    dotbracket_part = dotbracket[basepair[0]+1:basepair[1]]

    for element in dotbracket_part:
        if element == '(':
            count += 1
        if element == ')':
            count -= 1

    if count == 0:
        crossing = False

    log.debug('Crossing: baspair {}, count {}, crossing {}, {}'.format(basepair,count,crossing,dotbracket_part))

    return crossing


def is_taken(dotbracket, basepair):
    """
    Check if the current baspair is already taken (True|False)

    :param dotbracket: Current dotbracketline (list)
    :param baspair: Current basepair
    """
    taken = False

    open,close = basepair[0],basepair[1]
    if (dotbracket[open] != '.') or (dotbracket[close] != '.'):
        taken = True

    log.debug('Is taken: {} ({} {}, {} {})'.format(taken, open, dotbracket[open],close,dotbracket[close]))
    return taken


def create_db(basepairs, sequencelength,chainbreak):
    """
    Create a dotbracket string from the basepairlist

    :param basepairs: List of basepairs
    :param sequencelength: length of the DB
    """
    dotbracket = [['.']*sequencelength]
    dotbracket[0][chainbreak] = ' '
    dotbracketline = 0

    while len(basepairs) > 0:
        remaining_basepairs = []
        log.debug('NEW ROUND IN DB-CREATION')
        for position, basepair in enumerate(basepairs):
            log.debug('position {} basepair {}'.format(position, basepair))
            taken = is_taken(dotbracket[dotbracketline],basepair)
            crossing = is_crossing(dotbracket[dotbracketline],basepair)
            if (crossing == False) and (taken == False):
                log.debug('checkpoint passed {}'.format(basepair))
                dotbracket[dotbracketline][basepair[0]] = '('
                dotbracket[dotbracketline][basepair[1]] = ')'
            else:
                remaining_basepairs.append(basepair)
        basepairs = remaining_basepairs
        if len(basepairs) <= 0: break
        dotbracket.append(['.']*sequencelength)
        dotbracketline += 1
        dotbracket[dotbracketline][chainbreak] = ' '

    return dotbracket


def db2bps(db):
    """
    Get base pair list from db string.
    Possible as per DB-line there is no crosssing brackets allowed

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

def interfunction(db,chainbreak):
    """
    Isolate the interaction basepairs and measure the interaction length

    :param db: String of dotbracket
    :param chainbreak: (int) positon of the chainbreak between the two RNAs
    """
    #Find the interaction (use DBold)
    findinteractionline = dict()

    for nr, db_line in enumerate(db): #line includes only noncrossing bps
        log.debug("line: {}, db {}".format(nr, db_line))
        bp_count = 0
        interim_interaction_bps = db2bps(db_line)
        for pair in interim_interaction_bps: #counting interaction pairs
                if pair[0] < chainbreak and pair[1] > chainbreak:
                    bp_count += 1
        findinteractionline[nr] = [bp_count,interim_interaction_bps]
        log.debug(findinteractionline)

    interaction_list = list(sorted(findinteractionline.values()))[-1]
    interaction = interaction_list[0]= interaction_list[1]
    interaction.sort(key = lambda k: (k[0], -k[1]))

    log.debug('Interaction list, with interaction: {}'.format(interaction))

    return (interaction)


def main():
    parser = argparse.ArgumentParser(description='Prepare DB-files with expanding interaction')
    parser.add_argument ('-b', '--buffer', type= int, help='Buffer length',default=0)
    parser.add_argument ('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument ('-x', '--basepairlist', help='basepair list')
    parser.add_argument ('-d', '--dotbracket', help='dotbracket in SimRNA styles')
    parser.add_argument ('-n', '--nucleotides', help='related nucleotides')
    parser.add_argument ('-s', '--stepsize', type=int, help='e.g expanding of only 1/2 bps', default=1 )
    parser.add_argument ('-r', '--right', action='store_true', help='expand right')
    parser.add_argument ('-l', '--left', action='store_true', help='expand_left')
    parser.add_argument ('-o', '--output', help='name of the outputfile')
    parser.add_argument ('-t', '--target', help='end/target structure')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level =  logging.DEBUG)
        print('DEBUGGING MODE')
    else:
        logging.basicConfig(level = logging.WARNING)

    if (not args.right) and (not args.left): #default: both true
        right = True
        left = True
    else:
        right = args.right
        left = args.left

    log.debug('right {} left {}'.format(right, left))

    buffer = args.buffer
    Pairs = [['A', 'U'], ['C', 'G'], ['G', 'U'],['U', 'A'], ['G', 'C'], ['U', 'G']]

    #NUCLEOTIDES from file, SEQUENCELENGTH, CHAINBREAK
    nucleotides =  list()
    with open(args.nucleotides, 'r') as NCLFile:
        content = NCLFile.readlines()
        content = [s.rstrip() for s in content]
        for line in content:
            for index,element in enumerate(line):
                nucleotides.append(element)
                if not element.strip() or element == "&":
                    chainbreak = int(index)
        sequencelength = int(len(nucleotides))

    log.debug('ultimate ncl list with length {}, with chainbreak at {}'.format(sequencelength,chainbreak))
    log.debug(nucleotides)

    #BASEPAIRS
    basepairs = list()
    target_basepairs = list()
    if args.basepairlist is not None:
        with open(args.basepairlist, 'r') as BPFile:
            basepairs=json.load(BPFile)
        #basepairs.sort(key=itemgetter(0),reverse=False)
        basepairs.sort(key = lambda k: (k[0], -k[1]))
        log.debug('ultimate bp list: {}'.format(basepairs))

        #Starting Dotbracket structure
        dotbracket_old = create_db(basepairs,sequencelength,chainbreak)
        log.debug('old dotbracket')
        for line in dotbracket_old:
            log.debug(''.join(line))
    else:
        dotbracket_old=[]
        with open (args.dotbracket, 'r') as DBFile:
            for line in DBFile:
                dotbracket_old.append(line.rstrip('\n'))

        for db_line in dotbracket_old: #line includes only noncrossing bps
            interim_bp = (db2bps(db_line))
            for bp in interim_bp:
                basepairs.append(bp)
        basepairs.sort(key = lambda k: (k[0], -k[1]))
        log.debug('ultimate bp list: {}'.format(basepairs))
        log.debug('old dotbracket')
        for line in dotbracket_old:
            log.debug(''.join(line))

    interaction = interfunction(dotbracket_old,chainbreak)
    interactionlength= int(len(interaction))
    print('INTERACTIONLENGTH: {}'.format(interactionlength))

    #if there is a target structure given
    if args.target is not None:
        dotbracket_target=[]
        with open (args.target, 'r') as TARGETFile:
            for line in TARGETFile:
                dotbracket_target.append(line.rstrip('\n'))

        for target_line in dotbracket_target:
            interimtarget_bp = (db2bps(target_line))
            for bp in interimtarget_bp:
                target_basepairs.append(bp)

        target_basepairs.sort(key = lambda k: (k[0], -k[-1]))
        print('target bp list: {}'.format(target_basepairs))
        print('target dotbracket')
        for line in dotbracket_target:
            print(''.join(line))

        if basepairs == target_basepairs:
            print("Target structure reached")
            exit()

        target_interaction = interfunction(dotbracket_target,chainbreak)
        print(target_interaction)

    #new interaction + bufferzone
    start_left, end_left = list(), list()
    start_right, end_right = list(), list()

    target_i_left = list()
    target_i_right = list()

    if args.target is not None: #if a target structure is given
        interaction_firstpair = interaction[0][0]
        interaction_lastpair = interaction[-1][0]
        for nr, pair in enumerate(target_interaction):
            if pair[0] < interaction_firstpair:
                target_i_left.append(pair)
            if pair[0] > interaction_lastpair:
                target_i_right.append(pair)

        print('target interaction: left {} right {}'.format(target_i_left,target_i_right))
        if left and len(target_i_left) > 0:
            for i in range(target_i_left[-1][0]-buffer,interaction[0][0]):
                start_left.append(i)
            for i in range(interaction[0][-1]+1,target_i_left[-1][1]+buffer+1):
                end_left.append(i)

        if right and len(target_i_right) > 0:
            for i in range(interaction[-1][0]+1,target_i_right[0][0]+buffer+1):
                start_right.append(i)
            for i in range(target_i_right[0][1]-buffer,interaction[-1][1]):
                end_right.append(i)

        print('{}, {}, {}, {}'.format(start_left,end_left, start_right, end_right))

    else: #endless interaction expansion
        for i in range(buffer+1):
            i += 1 #one for zero and one for bracket + buffer
            if right and left:
                start_left.append(interaction[0][0]-i)
                end_left.append(interaction[0][1]+i)
                start_right.append(interaction[-1][0]+i)
                end_right.append(interaction[-1][1]-i)
            if right and not left:
                start_right.append(interaction[-1][0]+i)
                end_right.append(interaction[-1][1]-i)
            if left and not right:
                start_left.append(interaction[0][0]-i)
                end_left.append(interaction[0][1]+i)

            log.debug('{}, {}, {}, {}'.format(start_left,end_left, start_right, end_right))

    #Remove all brackets left
    if left:
        removing_basepairs = []
        start_left.sort()
        end_left.sort()
        for base in start_left:
            for pair in basepairs:
                if base in pair:
                    removing_basepairs.append(pair)
        for base in end_left:
            for pair in basepairs:
                if base in pair:
                    removing_basepairs.append(pair)
        for pair in removing_basepairs:
            basepairs.remove(pair)
        log.debug('New interaction positions left: start {} end {}'.format(start_left, end_left))

    if right:
        removing_basepairs = []
        start_right.sort()
        end_right.sort()
        for base in start_right:
            for pair in basepairs:
                if base in pair:
                    removing_basepairs.append(pair)
        for base in end_right:
            for pair in basepairs:
                if base in pair:
                    removing_basepairs.append(pair)
        for pair in removing_basepairs:
            basepairs.remove(pair)
        log.debug('New interaction positions right: start {} end {}'.format(start_right, end_right))



    if args.target is not None: #if a target structure is given
        if left and len(target_i_left) > 0:
            newbps=[target_i_left[-1][0],target_i_left[-1][1]]
            ncl = [nucleotides[newbps[0]],nucleotides[newbps[1]]]
            log.debug('new interaction left: {} {}, {} {}'.format(newbps[0],ncl[0],newbps[1],ncl[1]))
            #chainbraketest: chainbrake != in pairs
            if ncl in Pairs:
                basepairs.append(newbps)

        if right and len(target_i_right) > 0:
            newbps=[target_i_right[0][0],target_i_right[0][1]]
            ncl = [nucleotides[newbps[0]],nucleotides[newbps[1]]]
            log.debug('new interaction right: {} {}, {} {}'.format(newbps[0],ncl[0],newbps[1],ncl[1]))
            if ncl in Pairs:
                basepairs.append(newbps)

    else:
        #new interaction
        if left:
            newbps=[start_left[-1],end_left[0]]
            ncl = [nucleotides[start_left[-1]],nucleotides[end_left[0]]]
            log.debug('new interaction left: {} {}, {} {}'.format(newbps[0],ncl[0],newbps[1],ncl[1]))
            #chainbraketest: chainbrake != in pairs
            if ncl in Pairs:
                basepairs.append(newbps)

        if right:
            newbps=[start_right[0],end_right[-1]]
            ncl = [nucleotides[start_right[0]],nucleotides[end_right[-1]]]
            log.debug('new interaction right: {} {}, {} {}'.format(newbps[0],ncl[0],newbps[1],ncl[1]))
            if ncl in Pairs:
                basepairs.append(newbps)

    basepairs.sort(key = lambda k: (k[0], -k[1]))
    '''
        list = sort( [[5, 2], [2, 4], [4, 5]])
        list = [[2, 5], [2, 4], [4, 5]]
    '''
    log.debug('new basepairlists: {}'.format(basepairs))

    #WRITE THE NEW DOTBRACKET STRUCTURE
    dotbracket_new = create_db(basepairs,sequencelength,chainbreak)
    log.debug('new dotbracket')
    forprinting = list()
    for line in dotbracket_new:
        log.debug(''.join(line))
        fp = (''.join(line))
        forprinting.append(fp)

    #INTERACTION length
    interactionnew = interfunction(dotbracket_new,chainbreak)
    print('INTERACTIONLENGTH NEW: {}'.format(len(interactionnew)))

    #output only if there was an expansion
    if int(len(interactionnew)) > int(interactionlength):
        #WRITE OUTPUT
        with open(args.output, 'w') as out:
            out.write('\n'.join(forprinting))
            out.write('\n')
    else:
        print("No expansion possible any more!")
        exit()

if __name__ == "__main__":
    main()
