# RNA-Interaction-Workflow

Aufrufen mit (bevorzug Python3):
python expandinteraction.py -b 2 -x Interaction_Script/test4.bp -n Interaction_Script/test4.seq -i Interaction_Script/test4.in

Es geht die Interaktionsstelle einer RNA zu verlängern (wahlweise nach rechts, links oder beides (default))

BSP:
UUCGUACUCGCCAAAGUUGAA UCUUCUUCAACUUUGGCGAGUACGAAAAGA  --> RNA Sequence
((((.............)))) ((((.(((...............)))))))  --> Dotbracket Struktur
.....(((((((((((..... ..........))))))))))).........  --> Interaktion zwischen den zwei RNA Strängen
    RNA-STRANG 1                RNA-STRANG 2


Dotbracket Struktur: Punkt   = Nukleotid hat keine Bindung
                     Klammer = Nukleotid hat eine Bindung mit der jeweils anderen
                               e.g. 12    89       
                                    UUGUACAA
                                    ((....))
                               Basenpaare: [1,8] und [2,9]

-n ... File mit der Nukleotidsequence
-x ... Liste aller Basenpaare
-i ... alle Interactionspaare (zweite Dotbracketzeile)

Wenn nach rechts (R)/links (L) verlängert wird wird an folgenden Stellen ein neues Basenpaar eingefügt (vorausgesetzt sie sind komplementär (G-C, C-G, A-U, U-A)

((((.............)))) ((((.(((...............)))))))
....R(((((((((((L.... .........L)))))))))))R........


-b gibt den Buffer auf der rechten/linken seite von der Interaktion zur nächsten Klammer an (also wie viele Punkte bis zur nächsten Klammer) 
e.g. -b 2 (- sollen in der ersten Zeile dann . sein, also die neue Interaktion muss ein Punkt sein + die zwei Puffer)

(---...........---)) ((((.((---...........---))))))
....R(((((((((((L.... .........L)))))))))))R........


--> wenn man genau schaut, hat jetzt der erste Strang in der ersten Zeile eine Klammer zu viel ')'
Deswegen gibt es die Basenpaarliste - da werden immer die jeweile Basenpaare rausgelöscht, damit es dann am Ende so ausschaut:

(..................) ((((.((.................))))))
....(((((((((((((.... .........)))))))))))))........
