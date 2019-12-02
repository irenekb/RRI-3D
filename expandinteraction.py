#!usr/bin/env python
'''
python3 expandinteraction.py -b 2 -x Interaction_Script/test4.bp -n Interaction_Script/test4.seq -i Interaction_Script/test4.in

python expandinteraction.py -b 2 -x test4.bp -n test4.seq -i test4.in
'''

import logging
import argparse
#import csv
from collections import defaultdict
import json
from operator import itemgetter
import sys
import numpy as np
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


def main():
    parser = argparse.ArgumentParser(description='Prepare DB-files with expanding interaction')
    parser.add_argument ('-b', '--buffer', type= int, help='Buffer length',default=0)
    parser.add_argument ('-v', '--verbose', action='store_true', help='Be verbose')

    parser.add_argument ('-x', '--basepairlist', help='basepair list')
    parser.add_argument ('-n', '--nucleotides', help='related nucleotides')
    parser.add_argument ('-i', '--interaction', help='interacting basepairs')
    parser.add_argument ('-s', '--stepsize', type=int, help='e.g expanding of onlzy 1/2 bps', default=1 )
    parser.add_argument ('-r', '--right', action='store_true', help='expand right')
    parser.add_argument ('-l', '--left', action='store_true', help='expand_left')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level =  logging.DEBUG)
        print('DEBUGGING MODE')
    else:
        logging.basicConfig(level = logging.WARNING)

    if not args.right and not args.left:
        right = True
        left = True
    else:
        right = args.right
        left = args.left

    print('right {} left {}'.format(right, left))

    buffer = args.buffer
    Pairs = {"A":"U","U":"A","G":"C","C":"G"}
    #Pairs = [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')]

    #NUCLEOTIDES from file, SEQUENCELENGTH, CHAINBREAK
    nucleotides =  list()
    with open(args.nucleotides, 'r') as NCLFile:
        content = NCLFile.readlines()
        content = [s.rstrip() for s in content]
        for line in content:
            for index,element in enumerate(line):
                nucleotides.append(element)
                if not element.strip():
                    chainbreak = int(index)
        sequencelength = int(len(nucleotides))
    print('ultimate ncl list with length {}, with chainbrake at {}'.format(sequencelength,chainbreak))
    print(nucleotides)

    #BASEPAIRS
    basepairs = list()
    with open(args.basepairlist, 'r') as BPFile:
        basepairs=json.load(BPFile)
    #basepairs.sort(key=itemgetter(0),reverse=False)
    basepairs.sort(key = lambda k: (k[0], -k[1]))
    print('ultimate bp list: {}'.format(basepairs))


    #Starting Dotbracket structure
    dotbracket_old = create_db(basepairs,sequencelength,chainbreak)
    print('old dotbracket')
    for line in dotbracket_old:
        print(''.join(line))


    #Find the interaction (use DBold)
    findinteractionline = dict()
    for nr, db_line in enumerate(dotbracket_old):
        bp_count = 0
        interim_interaction_bps = db2bps(db_line)
        for pair in interim_interaction_bps:
                if pair[0] < chainbreak and pair[1] > chainbreak:
                    bp_count += 1
        findinteractionline[nr] = [bp_count,interim_interaction_bps]
        print (findinteractionline)

    interaction = max(findinteractionline.values()[1])
    #interaction.sort(key=itemgetter(0), reverse=False)
    interaction.sort(key = lambda k: (k[0], -k[1]))

    print('Interaction list, with interaction: {}'.format(interaction))

    #new interaction + bufferzone
    start_left, end_left = list(), list()
    start_right, end_right = list(), list()
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

    print(start_left,end_left, start_right, end_right)
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
        print('New interaction positions left: start {} end {}'.format(start_left, end_left,removing_basepairs))

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
        print('New interaction positions right: start {} end {}'.format(start_right, end_right,removing_basepairs))


    #new interaction
    if left:
        newbps=[start_left[-1],end_left[0]]
        ncl = [nucleotides[start_left[-1]],nucleotides[end_left[0]]]
        print('new interaction left: {} {}, {} {}'.format(newbps[0],ncl[0],newbps[1],ncl[1]))
        #chainbraketest: chainbrake != in pairs
        if ncl[0] in Pairs and ncl[1] == Pairs[ncl[0]]:
            basepairs.append(newbps)
    if right:
        newbps=[start_right[0],end_right[-1]]
        ncl = [nucleotides[start_right[0]],nucleotides[end_right[-1]]]
        print('new interaction right: {} {}, {} {}'.format(newbps[0],ncl[0],newbps[1],ncl[1]))
        if ncl[0] in Pairs and ncl[1] == Pairs[ncl[0]]:
        #if ncl in Pairs:
            basepairs.append(newbps)

    #basepairs = sorted(basepairs, key=lambda l: (len(l), l))
    basepairs.sort(key = lambda k: (k[0], -k[1]))
    '''
        list = sort( [[5, 2], [2, 4], [4, 5]])
        list = [[2, 5], [2, 4], [4, 5]]
    '''
    print('new basepairlists: {}'.format(basepairs))

    #WRITE THE NEW DOTBRACKET STRUCTURE
    dotbracket_new = create_db(basepairs,sequencelength,chainbreak)
    print('new dotbracket')
    for line in dotbracket_new:
        print(''.join(line))

if __name__ == "__main__":
    main()
