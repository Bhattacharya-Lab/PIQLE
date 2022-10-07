#Note: Run this program from the same directory where all the a3m files are

import os, sys, fnmatch, math
import numpy as np
import os,sys,shutil, re, subprocess, optparse, multiprocessing

parser = optparse.OptionParser()

parser.add_option('--tar', dest='targetName',
        default = '',    # default empty!
        help = 'Target name')
parser.add_option('--rec', dest='receptor',
        default = '',    # default empty!
        help = 'Receptor')
parser.add_option('--lig', dest='ligand',
        default = '',    # default empty!
        help = 'Ligand')
parser.add_option('--ch1', dest='chain1',
        default = '',    # default empty!
        help = 'Chain 1')
parser.add_option('--ch2', dest='chain2',
        default = '',    # default empty!
        help = 'Chain 2')
parser.add_option('--out', dest='outPath',
        default = '',    # default empty!
        help = 'Output path')

(options,args) = parser.parse_args()

targetName = options.targetName
rec = options.receptor
lig = options.ligand
chain1 = options.chain1
chain2 = options.chain2
outPath = options.outPath

PIQLE_path = 'change/to/your/current/directory'
toolDir = PIQLE_path + '/apps/glinter/external/'
ttFile = PIQLE_path + '/apps/glinter/tax_tree'

#generate nogap
print(toolDir + 'A3M_NoGap ' + rec + ' ' + outPath + '/' + rec + '.nogap')
os.system(toolDir + 'A3M_NoGap ' + rec + ' ' + outPath + '/' + rec + '.nogap')
os.system(toolDir + 'A3M_NoGap ' + lig + ' ' + outPath + '/' + lig + '.nogap')

#generate SpecBloc
os.system(toolDir + 'A3M_SpecBloc ' + outPath + '/' + rec + '.nogap ' +  ttFile + ' ' + outPath + '/' + rec + '.specbloc')
os.system(toolDir + 'A3M_SpecBloc ' + outPath + '/' + lig + '.nogap ' +  ttFile + ' ' + outPath + '/' + lig + '.specbloc')

#Concat MSA
os.system(toolDir + 'MSA_ConCat ' + outPath + '/' + rec + '.specbloc ' + outPath + '/' + lig + '.specbloc ' + outPath + '/' + targetName + '_' + chain1 + '_' + chain2 + '.concat.aln')

