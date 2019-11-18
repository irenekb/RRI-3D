#!usr/bin/env python
'''
python3 expandinteraction.py -b 2 -x Interaction_Script/test4.bp -n Interaction_Script/test4.seq -i Interaction_Script/test4.in

python expandinteraction.py -b 2 -x test4.bp -n test4.seq -i test4.in
'''

import logging
import argparse
from copy import deepcopy
#import csv
from collections import defaultdict
import json
from operator import itemgetter

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
        for position, basepair in enumerate(basepairs):
            log.debug(position, basepair)
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


def main():
    parser = argparse.ArgumentParser(description='Prepare DB-files with expanding interaction')
    parser.add_argument ('-b', '--buffer', type= int, help='Buffer length',default=0)
    parser.add_argument ('-v', '--verbose', action='store_true', help='Be verbose')

    parser.add_argument ('-x', '--basepairlist', help='basepair list')
    parser.add_argument ('-n', '--nucleotides', help='related nucleotides')
    parser.add_argument ('-i', '--interaction', help='interacting basepairs')
    parser.add_argument ('-s', '--stepsize', type=int, help='e.g expanding of onlzy 1/2 bps', default=1 )
    parser.add_argument ('-r', '--right', action='store_true', help='expand right',default=True)
    parser.add_argument ('-l', '--left', action='store_true', help='expand_left',default=True)

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level =  logging.DEBUG)
        print('DEBUGGING MODE')
    else:
        logging.basicConfig(level = logging.WARNING)

    buffer = args.buffer
    Pairs = {"A":"U","U":"A","G":"C","C":"G"}

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
    basepairs.sort(key=itemgetter(0),reverse=False)
    print('ultimate bp list: {}'.format(basepairs))


    #Starting Dotbracket structure
    dotbracket_old = create_db(basepairs,sequencelength,chainbreak)
    print('old dotbracket')
    for line in dotbracket_old:
        print(''.join(line))


    #INTERACTION BASEPAIRS
    interaction = list()
    with open(args.interaction, 'r') as INFile:
        interaction=json.load(INFile)
    interaction.sort(key=itemgetter(0), reverse=False)

    #new interaction + bufferzone
    start_left, end_left = list(), list()
    start_right, end_right = list(), list()
    for i in range(buffer+1):
        i += 1 #one for zero and one for bracket + buffer
        if args.right and args.left:
            start_left.append(interaction[0][0]-i)
            end_left.append(interaction[0][1]+i)
            start_right.append(interaction[-1][0]+i)
            end_right.append(interaction[-1][1]-i)
        if args.right and not args.left:
            start_right.append(interaction[-1][0]+i)
            end_right.append(interaction[-1][1]-i)
        if args.left and not args.right:
            start_left.append(interaction[0][0]-i)
            end_left.append(interaction[0][1]+i)

    if args.right:
        start_right.sort()
        end_right.sort()
    if args.left:
        start_left.sort()
        end_left.sort()

    print('Interaction list, with interaction: {}'.format(interaction))

    print('New interaction positions left: start {} end {}; right: start {} end {}'.format(start_left, end_left, start_right, end_right))

    #Remove all brackets left
    if args.left:
        for i in start_left:
            for k, v in enumerate(basepairs):
                if i in v:
                    basepairs.pop(k)
        for i in end_left:
            for k, v in enumerate(basepairs):
                if i in v:
                    basepairs.pop(k)
    #Remove all brackets right
    if args.right:
        for i in start_right:
            for k, v in enumerate(basepairs):
                if i in v:
                    basepairs.pop(k)
        for i in end_right:
            for k, v in enumerate(basepairs):
                if i in v:
                    basepairs.pop(k)

    print('cleaned basepairlists: {}'.format(basepairs))

    #new interaction
    if args.left:
        newbps=[start_left[-1],end_left[0]]
        ncl = [nucleotides[start_left[-1]],nucleotides[end_left[0]]]
        print('new interaction left: {} {}, {} {}'.format(newbps[0],ncl[0],newbps[1],ncl[1]))
        if ncl[0] in Pairs and ncl[1] == Pairs[ncl[0]]:
            basepairs.append(newbps)
    if args.right:
        newbps=[start_right[0],end_right[-1]]
        ncl = [nucleotides[start_right[0]],nucleotides[end_right[-1]]]
        print('new interaction right: {} {}, {} {}'.format(newbps[0],ncl[0],newbps[1],ncl[1]))
        if ncl[0] in Pairs and ncl[1] == Pairs[ncl[0]]:
            basepairs.append(newbps)
    #ToDo: action if chainbreak or it is not complementary

    basepairs = sorted(basepairs, key=lambda l: (len(l), l))
    '''
        list2 = [[4, 5, 2], [2, 5, 4], [2, 4, 5]]
        print(sorted(list2, key=lambda l: (len(l), l)))
        [[2, 4, 5], [2, 5, 4], [4, 5, 2]]
    '''
    print('new basepairlists: {}'.format(basepairs))

    #WRITE THE NEW DOTBRACKET STRUCTURE
    dotbracket_new = create_db(basepairs,sequencelength,chainbreak)
    print('new dotbracket')
    for line in dotbracket_new:
        print(''.join(line))

if __name__ == "__main__":
    main()
