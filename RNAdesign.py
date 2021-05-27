#!usr/bin/env python3
'''
RNAblueprint is a library to sample sequences that are compatible with multiple
structure constraints. This allows us to generate multi-stable RNAs, i.e. RNAs
that switch between multiple predefined structures.

This script allows to create designs for two interacting RNA structures.
The first structure describes the two seperate hairpins with a connection element (A)
the second structure should ensure the complementarity cleaveage of the two hairpins.
With the objective2 function every designed hearpin will be evaluated separately.

A selection of - default 5 - designs will be saved with the following filename:
    output + 'design'+ designnumber.seq


SAMPLE COMMAND:
python RNAdesign.py -i SimRNA_interaction/RNAdesign.test -o RNA

-i --input      secondary structure file in SimRNA format
-o --ouput      Ouputname
-n --number     number of designs
-s --selection  number of selected designs
'''

import RNAblueprint as rbp
import RNA
import math
import random
import sys
import pandas as pd
import argparse
import logging
log = logging.getLogger(__name__)

def testinput():
    '''
    Testinput if there is no input given.

    :return DOTBRACKET: secondary structure
    '''

    ################ TEST-INPUT ####################
    DOTBRACKET = ['(((((((.........))))))) (((((((.........)))))))',
                  '((((((((((((((((((((((( )))))))))))))))))))))))']

    return(DOTBRACKET)

def objective(sequence):
    '''
    Design bi-stable switches, ie.e sequences that will form either of the two
    structures with close to 50% probability. It uses RNAblueprint to sample
    sequences that can form both of these structures in principle.

    :param sequence
    '''

    # create fold compounds
    fc = RNA.fold_compound(sequence)
    # calculate ensemble free energy
    pf = fc.pf()[1]
    # energy of struct for 1st and 2nd structure
    eos1 = fc.eval_structure(STRUCTURES[1])
    eos0 = fc.eval_structure(STRUCTURES[0])
    # return objective function
    return (eos1 - pf) + (eos0 - pf);


def objective2(sequence,STRUCTURES,index):
    '''
    variant for 'interaction' between first and second half
    allows to spicify the stability in more detail between the two hairpins

    DIF. TO OBJECTIVE: optimize on both stems and ignore the second input structure
                        (only for complimentarity)

    :param sequence:
    :return: objective funtion
    '''
    # split the sequence/structure in the middle
    #chainbreak=int(len(sequence)/2)
    sequence1 = sequence[0:(index+5)]
    sequence2 = sequence[(index+5):]
    structure1 = STRUCTURES[0][0:(index+5)]
    structure2 = STRUCTURES[0][(index+5):]
    log.debug(STRUCTURES[0])
    log.debug(STRUCTURES[1])
    log.debug(sequence)
    log.debug(sequence1)
    log.debug(structure1)
    log.debug(sequence2)
    log.debug(structure2)

    # create the fold compounds of each part of the sequence
    fc1 = RNA.fold_compound(sequence1)
    fc2 = RNA.fold_compound(sequence2)
    # calculate the ensemble free energy of each part of the sequence
    pf1 = fc1.pf()[1]
    pf2 = fc2.pf()[1]
    # energy of each structural part
    eos1 = fc1.eval_structure(structure1)
    eos2 = fc2.eval_structure(structure2)

    return (eos1 - pf1) + (eos2 - pf2)


def main():
    '''
    Perform a simple optimization procedure using simulated annealing.
    The crucial part is the objective() function, which is designed such that it
    becomes minimal when the Boltzmann ensemble is dominated by the two target
    structures.
    '''

    parser = argparse.ArgumentParser(description='Design an RNA for a specific secondary structure')
    parser.add_argument ('-i', '--input', help='Secondary Structure - SimRNA Format')
    parser.add_argument ('-n', '--number', help='Number of Designs', default=10)
    parser.add_argument ('-s', '--selection', help='Number of selected Designs',
                            default=5)
    parser.add_argument ('-o', '--output',
                            help='Name of the outputfiles with consecutive numbers',
                            default="design")
    parser.add_argument ('-v', '--verbose', action='store_true', help='Be verbose')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level =  logging.DEBUG)
        print('DEBUGGING MODE')
    else:
        logging.basicConfig(level = logging.WARNING)

    NUBMER_DESIGNS = int(args.number)
    selection = int(args.selection)
    COOL_FACTOR = 0.98 #default
    STRUCTURES=[]

    if args.input:

        dotbracket = []
        with open (args.input, 'r') as INPUTFILE:
            for line in INPUTFILE:
                dotbracket.append(line.rstrip('\n'))
        log.debug(dotbracket)

    #if there is no input given, take the default one
    else:
        dotbracket = testinput()
        log.debug(dotbracket)

    #if there is a chainbreak connect the two chains
    index= dotbracket[0].index(" ")
    log.debug('chainbreak: {}'.format(index))
    for line in dotbracket:
        line = line[:index] + '..........' + line[index+1:]
        log.debug(line)
        STRUCTURES.append(line)

    #creat a sequence incl. the "chainbreaking-sequence"
    SEQUENCE_CONSTRAINTS= "N"* (len(dotbracket[0])-1) #-1 for the empty space
    SEQUENCE_CONSTRAINTS= SEQUENCE_CONSTRAINTS[:index] + 'AAAAAAAAAA' + SEQUENCE_CONSTRAINTS[index:]
    log.debug(SEQUENCE_CONSTRAINTS)
    erase_line = '\r\033[K'

    # print input
    print("\n".join(STRUCTURES) + "\n" + SEQUENCE_CONSTRAINTS + "\n")
    # generate dependency graph from input
    dg = rbp.DependencyGraphMT(STRUCTURES, SEQUENCE_CONSTRAINTS)

    sequence_df = pd.DataFrame(columns=['Sequence','Score'], index=None)

    for _ in range(0,NUBMER_DESIGNS):
        # reset and resample complete sequence
        dg.sample()
        # initialize best score,
        score = objective2(dg.get_sequence(),STRUCTURES,index);
        # initial temperature for cooling during simulated annealing
        temperature = 20

        # optimization loop with simulated annealing
        while temperature > 0.0001 :
            sys.stdout.write('%s%2.6f\t%2.6f\t%s' % (erase_line, temperature, score, dg.get_sequence()))
            sys.stdout.flush()
            temperature = temperature * COOL_FACTOR

            # sample a new sequence
            dg.sample_clocal()
            # calculate new score
            this_score = objective2(dg.get_sequence(),STRUCTURES,index)
            # evaluate probability with scores
            rand = random.uniform(0, 1)
            if (this_score-score) < 0:
                prob = 1
            else:
                prob = math.exp(-1*(this_score-score)/temperature)
            # compare and make decision
            if (rand <= prob):
                score = this_score
            else:
                # revert to previous sequence
                dg.revert_sequence()

        # finally collect and return the result
        sequence_df=sequence_df.append({'Sequence':dg.get_sequence(),'Score':score},ignore_index=True)
        sys.stdout.write('%s%s\t%2.6f\n' % (erase_line, dg.get_sequence(), score))
        sys.stdout.flush()

    #select the best files (selection)
    log.debug(sequence_df)
    sequence_df.sort_values('Score',ascending=True, inplace=True)
    sequence_df=sequence_df.head(selection)
    print(sequence_df)

    #save the best designs in several .seq file
    output = str(args.output)
    for n in range(0,selection):
        nstr= str(n+1)
        outputname = output + "design" + nstr + ".seq"

        #remove the "chainbreaking-sequence"
        sequence = str(sequence_df.iloc[n]['Sequence'])
        sequence= sequence[:index] +' '+ sequence[index+10:]
        with open(outputname, 'w') as out:
            out.write(sequence)

if __name__ == "__main__":
    main()
