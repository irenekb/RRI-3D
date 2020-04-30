Several RNA-RNA-interaction skripts


- expandinteraciton.py

Prepare DB-files with expanding interaction (python3 preferred)

Expand an interaction site between two RNA chains. 
The expansion can be at both sites of the interaction site (default), only the right or only the left interaction site.
Expansion of the interaction site only by complementary basepairs (G-C, C-G, A-U, U-A).

e.g.:
UUCGUACUCGCCAAAGUUGAA UCUUCUUCAACUUUGGCGAGUACGAAAAGA  --> RNA Sequence
((((.............)))) ((((.(((...............)))))))  --> Dotbracket structur (intra)
.....(((((((((((..... ..........))))))))))).........  --> Dotbracket structur (inter)
    RNA-CHAIN 1                RNA-CHAIN 2


Dotbracket Structure: Point   = free nucleotide
                      Bracket = basepair between two nucleotides
                               e.g. 12    89       
                                    UUGUACAA
                                    ((....))
                               basepair: [1,8] and [2,9]


> python expandinteraction.py -b 2 -x test4.bp -n test4.seq -i test4.in

-n ... file with the nucleotide sequence
-x ... list of all basepairs
-i ... all interacting basepairs (second line in the dotbracket file)

Expand the interaction right (-r) or left (-l):

((((.............)))) ((((.(((...............)))))))
....R(((((((((((L.... .........L)))))))))))R........


-b length of the buffer region, no intra- and interaction allowed, before and after the interaction site. 
e.g -b 2

(---............---)) ((((.((---...........---))))))
....R(((((((((((L.... .........L)))))))))))R........


(..................) ((((.((.................))))))
....(((((((((((((.... .........)))))))))))))........

-s ... how many nucleotides should be added to the interaction (on one site)



- SSalignment.py

Compare all secondary structures files (calculated via SimRNA – SimRNA-style) with the start ss-sequence, the constrained ss-sequence and with each other.
For comparing the energy the SimRNA trafl file is used. 

The first output is a csv-file with the following information regarding one specific random number seed of a full run:
number, sequence, count_constraint, count_start, count_before, constancy, dif_constraint, dif_start, dif_before, bp, time, energy_values_plus_restraint_score, energy_value, current_temp, interaction, len_interaction, count_interaction_constraint, dif_interaction_cc

The second output is a csv.file with all unique structures over several (offered) random number seeds of a full run:
sequence, count_how_often, count_constraint, count_start, dif_constraint, dif_start, bp, bpstr, interaction, len_interaction, count_interaction_constraint, dif_interaction_cc

-m ... overwrite ('w') or append ('a') a random number seed run to the second output. 

> python3 SSalignment.py p /place/with/all/ss-quences -i ss-constrain -c ssstart -o firstoutput.csv -u secondoutput.csv -m 'w' -t traflfile



- continoussearch.py

Parse over all SSalignment files (individual runs and the overview).
Looking for the most common secondary  structure in the overview file.
Look for this secondary structure in all individual runs and separate them (max_file).
The structure with the best energy (3D) is the one for the next constrained run.


instead of the most common secondary structure: 
--force: Find the secondary structure that fits the given constrain best

--interaction: Find the secondary structure with an interaction most similar to the contrained one

general options:
- p       ... path to input files
-– print  ... print a csv-file with all minEnergy relevant files
–- first  ... verify the first line
–- second ... verify the second line

> python3 continoussearch.py -p path/to/inputfiles --print --first 'CopStems_00.ss' --second 'CopStems_00_00_000000.ss' -f


- To get the pdb for the next run:
SimRNA_trafl2pdbsstructure.pdb trajectory.trafl {list} AA
AA   ... all atom reconstruction
list ... frame in the list

