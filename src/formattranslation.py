#!usr/bin/env python3
'''
Create a fasta-file format with
    - several RNAdesigns
    - secondary structrue in simRNA format
    - structure name

TESTINPUT
python formattranslation.py -p PATHtoINPUTFILES -n NAMEofINPUTFILES -c 100
'''

import logging
import argparse

log = logging.getLogger(__name__)

def db2bps(db):
    """
    Get base pair list from db string.

    :param db: String of dotbracket
    """
    opening_bps = []
    bps = []
    chainbreak = int()

    for i in range(len(db)):
        if db[i] in ['(','[','{','<']:
            opening_bps.append(i)
        elif db[i] in [')',']','}','>']:
            bps.append([opening_bps.pop(),i]) # list of lists
        elif not db[i].strip() or db[i] =='&':
            chainbreak = i
        elif db[i] =='.':
            continue
        else:
            log.warning('Wrong character in dotbracketstring: {}'.format(db[i]))

    return bps, chainbreak

def interfunction(db):
    """
    Isolate the interaction basepairs and measure the interaction length

    :param db: String of dotbracket
    :param chainbreak: (int) positon of the chainbreak between the two RNAs
    """
    findinteractionline = dict()

    for nr, db_line in enumerate(db): #line includes only noncrossing bps
        bp_count = 0
        interim_interaction_bps,chainbreak = db2bps(db_line)
        for pair in interim_interaction_bps: #counting interaction pairs
            if pair[0] < chainbreak and pair[1] > chainbreak:
                bp_count += 1
        findinteractionline[nr] = [bp_count,interim_interaction_bps]

    interaction_list = list(sorted(findinteractionline.values()))[-1]
    interaction = interaction_list[0]= interaction_list[1]
    interaction.sort(key = lambda k: (k[0], -k[1]))

    return (interaction)


def transformernwin(bps,interaction,length,chainbreak):
    """
    Transform both bps lists (intra/intermolecular bps) into a fasta dotbracket format

    :param bps: list of basepairs
    :param interaction: list of interacting basepairs
    """
    dotbracket_ernwin=['.']*(length-1)

    for bp in bps:
        dotbracket_ernwin[bp[0]]='('
        dotbracket_ernwin[bp[1]]=')'

    for bp in interaction:
        dotbracket_ernwin[bp[0]]='['
        dotbracket_ernwin[bp[1]]=']'

    dotbracket_ernwin[chainbreak]='&'

    return dotbracket_ernwin


def main():
    parser = argparse.ArgumentParser(description='Prepare DB-files with expanding interaction')
    parser.add_argument ('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument ('-n', '--name', type=str, help='base name')
    parser.add_argument ('-c', '--count', type=int, help='number of samples')
    parser.add_argument ('-p', '--path', help='Path to Inputfiles')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level =  logging.DEBUG)
        print('DEBUGGING MODE')
    else:
        logging.basicConfig(level = logging.WARNING)

    path = args.path
    name = args.name
    sstructure =path+"/"+name+'_0.ss'
    log.debug(sstructure)

    dotbracket_simrna=[]
    dotbracket_ernwin=[]
    bp_simrrna=[]

    with open(sstructure, 'r') as SSTRUCTURE:
        for line in SSTRUCTURE:
            dotbracket_simrna.append(line.rstrip('\n'))
            length = int(len(line))

    log.debug('sequencelength: {}'.format(length))
    log.debug('simrna dotbracket')
    for line in dotbracket_simrna:
        log.debug(line)

    for db_line in dotbracket_simrna:
        interim_bps,chainbreak = (db2bps(db_line))
        for bp in interim_bps:
             bp_simrrna.append(bp)
        bp_simrrna.sort(key = lambda k: (k[0], -k[1]))

    log.debug('ultimate bp list: {}'.format(bp_simrrna))

    #interaction = interfunction(bp_simrrna,chainbreak)
    interaction = interfunction(dotbracket_simrna)
    log.debug('interaction {}'.format(interaction))
    log.debug('interactionlength: {}'.format(len(interaction)))

    dotbracket_ernwin=transformernwin(bp_simrrna,interaction,length,chainbreak)

    for n in range(args.count):
        n=str(n+1)
        SEQUENCEFILE = path+"/"+name+'design'+n+'.seq'
        OUTPUTNAME = path+"/"+name+'design'+n+'.fa'
        log.debug('Sequencefile: {} Outputname: {}'.format(SEQUENCEFILE,OUTPUTNAME))
        with open (SEQUENCEFILE, 'r') as FILE:
            for line in FILE:
                sequence=str(line).replace(' ', '&') #for fasta output

        output = open(OUTPUTNAME, 'w')
        output.writelines('>'+name+n+'\n')
        output.writelines(sequence+'\n')
        output.writelines(dotbracket_ernwin)


if __name__ == "__main__":
    main()
