#!usr/bin/env python3
'''
RNAblueprint is a library to sample sequences that are compatible with multiple
structure constraints. This allows us to generate multi-stable RNAs, i.e. RNAs
that switch between multiple predefined structures.
The first structure describes the two seperate hairpins with a connection element (A)
the second structure should ensure the complementarity cleaveage of the two hairpins.
'''


# without these two lines "import RNAblueprint" will fail
from ctypes import *
lib1 = cdll.LoadLibrary('/usr/local/lib/libRNAblueprint.so')

import RNAblueprint as rbp
import RNA
import math
import random
import sys

################ INPUT ####################
STRUCTURES =          ['(((((((.........)))))))..........(((((((.........)))))))',
                       '(((((((((((((((((((((((..........)))))))))))))))))))))))']
SEQUENCE_CONSTRAINTS = 'NNNNNNNNNNNNNNNNNNNNNNNAAAAAAAAAANNNNNNNNNNNNNNNNNNNNNNN'
NUBMER_DESIGNS = 50
COOL_FACTOR = 0.98

###########################################
erase_line = '\r\033[K'


def objective(sequence):
    '''
    Design bi-stable switches, ie.e sequences that will form either of the two
    structures with close to 50% probability. It uses RNAblueprint to sample
    sequences that can form both of these structures in principle.

    :param sequence:
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


def objective2(sequence):
    '''
    variant for 'interaction' between first and second half
    allows to spicify the stability in more detail between the two hairpins

    DIF. TO OBJECTIVE: optimize on both stems and ignore the second input structure
                        (only for complimentarity)

    :param sequence:
    :return: objective funtion
    '''
    # split the sequence/structure in the middle
    sequence1 = sequence[0:len(sequence)/2]
    sequence2 = sequence[len(sequence)/2:]
    structure1 = STRUCTURES[0][0:len(sequence)/2]
    structure2 = STRUCTURES[0][len(sequence)/2:0]

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

    # print input
    print("\n".join(STRUCTURES) + "\n" + SEQUENCE_CONSTRAINTS + "\n")
    # generate dependency graph from input
    dg = rbp.DependencyGraphMT(STRUCTURES, SEQUENCE_CONSTRAINTS)

    for _ in range(0,NUBMER_DESIGNS):
        # reset and resample complete sequence
        dg.sample()
        # initialize best score
        score = objective(dg.get_sequence());
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
            this_score = objective(dg.get_sequence())
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

        # finally return the result
        sys.stdout.write('%s%s\t%2.6f\n' % (erase_line, dg.get_sequence(), score))
        sys.stdout.flush()

if __name__ == "__main__":
    main()
