#!usr/bin/env python

import logging
import argparse
import copy
import csv
from collections import defaultdict
import pandas

import forgi
#import forgi.graph.bulge_graph as fgb
import forgi.utilities.commandline_utils as fuc
#import forgi.threedee.utilities.vector as ftuv
#from forgi.utilities.exceptions import GraphIntegrityError
#from forgi.utilities.exceptions import GraphConstructionError
#from numbered_dotbracket import NumberedDotbracket
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.vector as ftuv
import itertools as it

import sys
from optparse import OptionParser


log = logging.getLogger(__name__)

cg = forgi.load_rna('/home/irene/Documents/Studium/PhD/m_virus/BVDV/firstSimRNA/IRES5_180-204_357-425_final_A/simulation_01/step003000_BESTTRANSLATION.pdb')
print (list(cg.iloop_iterator()))
for iloop in cg.iloop_iterator():
     conn = cg.connections(iloop)
     # conn contains two values ['s0', 's1']
     angle = ftuv.vec_angle(cg.coords.get_direction(conn[0]), cg.coords.get_direction(conn[1]))
     print(iloop, angle)
