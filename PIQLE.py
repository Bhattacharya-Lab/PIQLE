#!/usr/bin/python
import os,sys,shutil, math
import optparse, argparse    # for option sorting
import csv
from decimal import *
from multiprocessing import Process
import multiprocessing
from threading import Thread
import subprocess
from itertools import combinations
import numpy as np
import time

import torch as th
import os, sys, fnmatch
import os,sys,shutil, re, subprocess, optparse, multiprocessing
import dgl
import dgl.function as fn
import torch
import torch.nn as nn
import torch.nn.functional as F


#------------------------configure------------------------------#
#Checking configuration                                         #
#---------------------------------------------------------------#
configured = 0
PIQLE_path = 'change/to/your/current/directory'
if(configured == 0 or not os.path.exists(PIQLE_path + '/apps/dssp') or
           not os.path.exists(PIQLE_path + '/apps/calNf_ly') or
           not os.path.exists(PIQLE_path + '/apps/stride') or
           not os.path.exists(PIQLE_path + '/scripts/distance_generator.py') or
           not os.path.exists(PIQLE_path + '/scripts/msa_concat.py') or
           not os.path.exists(PIQLE_path + '/scripts/msa_processor.py') or
           not os.path.exists(PIQLE_path + '/scripts/orientation_generator.py')):
	print("\nError: not yet configured!\nPlease configure as follows\n$ cd PIQLE\n$ python config.py\n")
	exit(1)

dssp_path = PIQLE_path + '/apps/dssp'
stride_path = PIQLE_path + '/apps/stride'
neff_generator = PIQLE_path + '/apps/calNf_ly'
orientation_path = PIQLE_path + '/scripts/orientation_generator.py'
distance_path = PIQLE_path + '/scripts/distance_generator.py'
msa_concat = PIQLE_path + '/scripts/msa_concat.py'
process_msa = PIQLE_path + '/scripts/msa_processor.py'

parser = argparse.ArgumentParser()

parser._optionals.title = "Arguments"
# take input arguments
parser.add_argument('--tgt', dest='targetName',
        default = '',    # default empty!
        help = 'Target name')

parser.add_argument('--seq', dest='fastaFile',
        default = '',    # default empty!
        help = 'Fasta file')

parser.add_argument('--dec', dest='decoyDir',
        default = '',    # default empty!
        help = 'Complex decoy directory')

parser.add_argument('--ch', dest='chainFile',
	default = '',    # default empty!
	help = 'Chain file')

parser.add_argument('--msa1', dest='inMsa1',
	default = '',    # default empty!
	help = 'MSA1: Multiple Sequence Alignment of chain 1')

parser.add_argument('--msa2', dest='inMsa2',
	default = '',    # default empty!
	help = 'MSA2: Multiple Sequence Alignment of chain 2')

parser.add_argument('--a3m1', dest='inA3m1',
	default = '',    # default empty!
	help = 'A3M of chain1')

parser.add_argument('--a3m2', dest='inA3m2',
	default = '',    # default empty!
	help = 'A3M of chain2')

parser.add_argument('--out', dest='outdir',
	default = '',    # default empty!
	help = 'Output directory.')

if len(sys.argv) < 5:
        parser.print_help(sys.stderr)
        sys.exit(1)

options = parser.parse_args()

targetName = options.targetName
fasta = options.fastaFile
decoyPath = options.decoyDir
chainFile = options.chainFile
msa1 = options.inMsa1
msa2 = options.inMsa2
a3m1 = options.inA3m1
a3m2 = options.inA3m2
outPath = options.outdir

decoys = os.listdir(decoyPath)
working_path = os.getcwd()
intFDist = 10

#-----create necessary directories-------#
#                                        #
#----------------------------------------#

#root output directory
if not(os.path.exists(outPath)):
    os.system('mkdir -p ' + outPath)
    
#directory for unbound decoys
if not(os.path.exists(outPath + '/decoys/unbound/')):
    os.system('mkdir -p ' + outPath + '/decoys/unbound/')
    
#directory for fasta files
if not(os.path.exists(outPath + '/fasta/unbound/')):
    os.system('mkdir -p ' + outPath + '/fasta/unbound/')

#------------read chain ids--------------#
#                                        #
#----------------------------------------#
chainIds = []
with open(chainFile) as cFile:
    for line in cFile:
        tmp = line.split()
        for i in range(len(tmp)):
            chainIds.append(tmp[i].strip())
        break

#-----------split fasta files------------#
#                                        #
#----------------------------------------#
chainNo = 0
with open(fasta) as fFile:
    for line in fFile:
        if(line[0] == '>'):
            outFasta = open(outPath + '/fasta/unbound/' + targetName + '_' + chainIds[chainNo] + '.fasta', 'w')
            chainNo += 1
        outFasta.write(line)
    outFasta.close()

#-----------split decoy files------------#
#                                        #
#----------------------------------------#
for d in range(len(decoys)):
    #for each pdb file
    lines = []
    with open(decoyPath + '/' + decoys[d]) as pFile:
        for line in pFile:
            if(len(line.split())>0 and line.split()[0] == "ATOM"):
                lines.append(line)

    #split the pdb for each chain and save to unbound directory
    for c in range(len(chainIds)):
        if not(os.path.exists(outPath + '/decoys/unbound/' + targetName + '_' + chainIds[c])):
            os.system('mkdir -p ' + outPath + '/decoys/unbound/' + targetName + '_' + chainIds[c])
        outFile = open(outPath + '/decoys/unbound/' + targetName + '_' + chainIds[c] + '/' + decoys[d].split('.pdb')[0] + '_' + chainIds[c] + '.pdb', 'w')
        for l in range(len(lines)):
            if(lines[l][21:(21+1)].strip() == chainIds[c]):
                outFile.write(lines[l])
        outFile.close()

def concatMSA():
    os.chdir(working_path)
    if not(os.path.exists(outPath + '/msa')):
        os.system('mkdir ' + outPath + '/msa')
    os.chdir(outPath + '/msa/')

    for c in range(len(chainIds)):
        os.system('cp ' + msa1 + ' ./')
        os.system('cp ' + msa2 + ' ./')
        os.system('cp ' + a3m1 + ' ./')
        os.system('cp ' + a3m2 + ' ./')
    
    for comp in combinations(chainIds, 2):
        chains = list(comp)
        recTemp = os.path.basename(os.path.normpath(a3m1))
        ligTemp = os.path.basename(os.path.normpath(a3m2))
        print('python ' + msa_concat + ' --tar ' + targetName + ' --rec ' + a3m1 + ' --lig ' +
                  a3m2 + ' --ch1 ' + chainIds[0] + ' --ch2 ' + chainIds[1] + ' --out ./')
        os.system('python ' + msa_concat + ' --tar ' + targetName + ' --rec ' + recTemp + ' --lig ' +
                  ligTemp + ' --ch1 ' + chainIds[0] + ' --ch2 ' + chainIds[1] + ' --out ./')
        
    os.chdir(working_path)

def calculateNeff():
    os.chdir(working_path)
    total = (len(chainIds) * (len(chainIds) - 1)) / 2

    con_msa = os.listdir(outPath + '/msa/')

    while(1):
        done = 0
        for con in range(len(con_msa)):
            if((outPath + '/msa/' + con_msa[con]).endswith('.concat.aln') and os.path.getsize(outPath + '/msa/' + con_msa[con]) > 0):
                done += 1
        if(done == total):
            break

    if not(os.path.exists(outPath + '/neff/bound')):
        os.system('mkdir -p ' + outPath + '/neff/bound')

    if not(os.path.exists(outPath + '/neff/unbound')):
        os.system('mkdir -p ' + outPath + '/neff/unbound')

    #process bound neff
    os.system('python ' + process_msa + ' --aln ' + outPath + '/msa/ --out ' + outPath + '/neff/bound')
    os.chdir(outPath + '/neff/bound')
    alns = os.listdir(outPath + '/neff/bound/')
    for a in range(len(alns)):
        if(alns[a].endswith('.concat.aln')):
            os.system(neff_generator + ' ' + outPath + '/neff/bound/' + alns[a] + ' 0.8 > ' + alns[a].split('.')[0] + '.neff')

    os.chdir(working_path)
    #process unbound neff
    os.chdir(outPath + '/neff/unbound')
    alns = os.listdir(outPath + '/msa/')
    for a in range(len(alns)):
        if(alns[a].endswith('.aln') and 'concat' not in alns[a]):
            print(neff_generator + ' ' + outPath + '/msa/' + alns[a] + ' 0.8 > ' + alns[a].split('.')[0] + '.neff')
            os.system(neff_generator + ' ' + outPath + '/msa/' + alns[a] + ' 0.8 > ' + alns[a].split('.')[0] + '.neff')

    os.chdir(working_path)

    
def runDSSP():
    os.chdir(working_path) #just in case
    for c in range(len(chainIds)):
        if not(os.path.exists(outPath + '/dssp/unbound/' + targetName + '_' + chainIds[c])):
            os.system('mkdir -p ' + outPath + '/dssp/unbound/' + targetName + '_' + chainIds[c])
        os.chdir(outPath + '/dssp/unbound/' + targetName + '_' + chainIds[c])
        decoys = os.listdir(outPath + '/decoys/unbound/' + targetName + '_' + chainIds[c])
        for d in range(len(decoys)):
            dssp_ret_code = os.system(dssp_path + ' -i ' + outPath + '/decoys/unbound/' + targetName + '_' + chainIds[c] + '/' + decoys[d] +
                      ' -o ' + decoys[d].split('.pdb')[0] + '.dssp')

            if(dssp_ret_code != 0):
                print("DSSP failed to run. Running STRIDE for " + outPath + '/decoys/unbound/' + targetName + '_' + chainIds[c] + '/' + decoys[d])
                os.system(stride_path +" " + outPath + '/decoys/unbound/' + targetName + '_' + chainIds[c] + '/' + decoys[d] + ">" +
                        decoys[d].split('.pdb')[0] + '.stride')
        os.chdir(working_path)
            

def generatePairs():
    #it will run only for targets with more than 2 chains
    if(len(chainIds) > 1):
        #make complex PDBs with 2 chains
        for comp in combinations(chainIds, 2):
            chains = list(comp)
            if not(os.path.exists(outPath + '/decoys/bound/' + targetName + '_' + chains[0] + '_' + chains[1])):
                os.system('mkdir -p ' + outPath + '/decoys/bound/' + targetName + '_' + chains[0] + '_' + chains[1])

            #for each decoys
            for d in range(len(decoys)):
                lines = []
                #for each chain
                chainNo = 0
                for c in chains:
                    with open(decoyPath + '/' + decoys[d]) as dFile:
                        for line in dFile:
                            if(len(line.split())>0 and line.split()[0] == "ATOM"):
                                if(line[21:(21+1)].strip() == chains[chainNo]):
                                    lines.append(line)
                    chainNo += 1

                outFile = open(outPath + '/decoys/bound/' + targetName + '_' + chains[0] + '_' + chains[1] + '/' + decoys[d], 'w')
                for l in range(len(lines)):
                    outFile.write(lines[l])
                outFile.close()

def generateOrientation():
    boundDecDir = os.listdir(outPath + '/decoys/bound/')
    for b in range(len(boundDecDir)):
        if(os.path.isdir(outPath + '/decoys/bound/' + boundDecDir[b])):
            if not(os.path.exists(outPath + '/orientation/' + boundDecDir[b])):
                os.system('mkdir -p ' + outPath + '/orientation/' + boundDecDir[b])

            chain1 = boundDecDir[b].split('_')[1]
            chain2 = boundDecDir[b].split('_')[2]
            os.system('python ' + orientation_path + ' -d ' + outPath + '/decoys/bound/' + boundDecDir[b] + ' -l ' + chain1 + ' -r ' +
                      chain2 + ' -o ' + outPath + '/orientation/' + boundDecDir[b])

def generateDistance():
    boundDecDir = os.listdir(outPath + '/decoys/bound/')
    for b in range(len(boundDecDir)):
        if(os.path.isdir(outPath + '/decoys/bound/' + boundDecDir[b])):
            if not(os.path.exists(outPath + '/distance/' + boundDecDir[b])):
                os.system('mkdir -p ' + outPath + '/distance/' + boundDecDir[b])

            chain1 = boundDecDir[b].split('_')[1]
            chain2 = boundDecDir[b].split('_')[2]

            os.system('python ' + distance_path + ' -d ' + outPath + '/decoys/bound/' + boundDecDir[b] + ' -l ' + chain1 + ' -r ' +
                      chain2 + ' -o ' + outPath + '/distance/' + boundDecDir[b])

#----------feature generation------------#
#                                        #
#----------------------------------------#
def get8to3ss(ss_parm):
    eTo3=""
    if (ss_parm == "H" or ss_parm == "G" or ss_parm == "I"):
            eTo3="H"
    elif(ss_parm == "E" or ss_parm == "B"):
            eTo3="E"
    else:
            eTo3="C"
    return eTo3

def get8to3ssOneHot(ss_parm):
    eTo3=""
    if (ss_parm == "H" or ss_parm == "G" or ss_parm == "I"):
            eTo3="1 0 0"
    elif(ss_parm == "E" or ss_parm == "B"):
            eTo3="0 1 0"
    else:
            eTo3="0 0 1"
    return eTo3

def edgeFeatEncoding(dist):
    dist = float(dist) #just in case
    feat = ''
    if(dist <= 2):
        feat = '1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
    if(dist > 2 and dist <= 2.5):
        feat = '0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
    if(dist > 2.5 and dist <= 3):
        feat = '0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
    if(dist > 3 and dist <= 3.5):
        feat = '0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0'
    if(dist > 3.5 and dist <= 4):
        feat = '0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0'
    if(dist > 4 and dist <= 4.5):
        feat = '0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0'
    if(dist > 4.5 and dist <= 5):
        feat = '0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0'
    if(dist > 5 and dist <= 5.5):
        feat = '0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0'
    if(dist > 5.5 and dist <= 6):
        feat = '0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0'
    if(dist > 6 and dist <= 6.5):
        feat = '0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0'
    if(dist > 6.5 and dist <= 7):
        feat = '0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0'
    if(dist > 7 and dist <= 7.5):
        feat = '0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0'
    if(dist > 7.5 and dist <= 8):
        feat = '0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0'
    if(dist > 8 and dist <= 8.5):
        feat = '0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0'
    if(dist > 8.5 and dist <= 9):
        feat = '0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0'
    if(dist > 9 and dist <= 9.5):
        feat = '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0'
    if(dist > 9.5 and dist <= 10):
        feat = '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1'
    return feat

def getSolAccOneHot(amnAcidParam, saValParam):
    saNorm="";
    aaSA=[];
    aaSaVal=[];
    aaSA=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'];
    aaSaVal=[115, 135, 150, 190, 210, 75, 195, 175, 200, 170, 185, 160, 145, 180, 225, 115, 140, 155, 255, 230];
    k=0
    while k<len(aaSA):
        if(amnAcidParam==aaSA[k]):
            sVal=((100 * float(saValParam)) / aaSaVal[k]);
            if(sVal<25):
                saNorm="1 0"
            elif(sVal>=25):
                saNorm="0 1"
            break;
        else:
            k+=1;
    return saNorm

def aaGroup(aa):
    aa = aa.capitalize().strip()
    if(aa == 'A' or aa == 'V' or aa == 'L' or aa == 'I' or aa == 'P' or aa == 'F' or aa == 'M' or aa == 'W'): #non-polar
        label = '1 0 0 0 0'
    elif(aa == 'G' or aa == 'S' or aa == 'T' or aa == 'C' or aa == 'Y' or aa == 'N' or aa == 'Q'): #polar
        label = '0 1 0 0 0'
    elif(aa == 'D' or aa == 'E'): #acidic
        label = '0 0 1 0 0'
    elif(aa == 'K' or aa == 'R' or aa == 'H'): #basic
        label = '0 0 0 1 0'
    else:
        label = '0 0 0 0 1'

    return label

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def sigmoid(x):
      if x < 0:
            return 1 - 1/(1 + math.exp(x))
      else:
            return 1/(1 + math.exp(-x))

def contains_number(str):
    return any(char.isdigit() for char in str)

def get_unique_list(in_list):
        if isinstance(in_list,list):
                return list(set(in_list))

def get3to1aa(aa):
    dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    if(aa not in dict):
        return 'X'
    else:
        return dict[aa]

def get_neff(neffFile):
    neff = 0
    with open(neffFile) as n_file:
        for line in n_file:
            tmp = line.split()
            if(len(tmp) > 0):
                x=np.array(tmp)
                x=np.asfarray(x, float)
                neff = sum(x) / len(x)

    return neff

def getSequence(pdb):
    orderedSeqNum = []
    originalSeqNum = []
    residueList = []
    start = 1
    
    tmp_residue_list = []
    start_end_ResNo = []
    prev_res_no = -1
    with open(pdb) as file:                                                   
        for line in file:                                                                                 
            if(line[0:(0+4)]=="ATOM"):
                if(prev_res_no != int(line[22:(22+4)].strip())):
                    residueList.append(get3to1aa(line[17:(17+3)].strip()))
                    originalSeqNum.append(int(line[22:(22+4)].strip()))
                    orderedSeqNum.append(start)
                    start += 1
                    
                prev_res_no = int(line[22:(22+4)].strip())
                
    return orderedSeqNum, originalSeqNum, residueList

def fastaLength(fastaFile):
    length  = 0
    with open(fastaFile) as fFile:
        for line in fFile:
            if(line[0] == '>'):
                continue
            length += len(line.strip())
    return length

def generateFeatures():
    total = (len(chainIds) * (len(chainIds) - 1)) / 2
    while(1):
        done_unbound = 0
        for c in range(len(chainIds)):
            if(os.path.exists(outPath + '/neff/unbound/' + targetName + '_' + chainIds[c] + '.neff') and
               os.path.getsize(outPath + '/neff/unbound/' + targetName + '_' + chainIds[c] + '.neff') > 0):
                done_unbound += 1

        done_bound = 0
        for comp in combinations(chainIds, 2):
            chains = list(comp)
            if(os.path.exists(outPath + '/neff/bound/' + targetName + '_' + chains[0] + '_' + chains[1] + '.neff') and
               os.path.getsize(outPath + '/neff/bound/' + targetName + '_' + chains[0] + '_' + chains[1] + '.neff') > 0):
                done_bound += 1
                
        if(done_unbound == len(chainIds) and done_bound == total):
            break
    print("Generating features...")
    #for each dimer
    complexDir = os.listdir(outPath + '/distance')
    #for each distance file
    for c in range(len(complexDir)):
        if not(os.path.exists(outPath + '/features/' + complexDir[c])):
            os.system('mkdir -p ' + outPath + '/features/' + complexDir[c])

        print(complexDir[c])
        files_rr = os.listdir(outPath + '/distance/' + complexDir[c])
        for f in range(len(files_rr)):
            if((outPath + '/distance/' + complexDir[c] + '/' + files_rr[f]).endswith('.rr')):
                print(outPath + '/distance/' + complexDir[c] + '/' + files_rr[f])
                decLines = []
                with open(outPath + '/distance/' + complexDir[c] + '/' + files_rr[f]) as dFile:
                    for line in dFile:
                        if(float(line.split()[2]) < intFDist):
                            decLines.append(line)

            #---get lines from dssp or stride---#
            #                                   #
            #-----------------------------------#
            dsspChain1Found = 0
            dsspChain2Found = 0

            chain1Dssp = []
            line1Count = 0
            
            chain2Dssp = []
            line2Count = 0

            chain1Stride = []
            chain2Stride = []

            if(os.path.exists(outPath + '/dssp/unbound/' + targetName + '_' + complexDir[c].split('_')[1] + '/' +
                              files_rr[f].split('.rr')[0] + '_' + complexDir[c].split('_')[1] + '.dssp')):

                #read dssp lines

                with open(outPath + '/dssp/unbound/' + targetName + '_' + complexDir[c].split('_')[1] + '/' +
                              files_rr[f].split('.rr')[0] + '_' + complexDir[c].split('_')[1] + '.dssp') as dssp1File:
                    for line in dssp1File:
                        if (line1Count<1):
                            if (line[2:(2+1)] == '#'):
                                line1Count += 1
                                continue
                        if(line1Count > 0):
                            if(len(line) > 0):
                                chain1Dssp.append(line)
                    dsspChain1Found = 1

            elif(os.path.exists(outPath + '/dssp/unbound/' + targetName + '_' + complexDir[c].split('_')[1] + '/' +
                              files_rr[f].split('.rr')[0] + '_' + complexDir[c].split('_')[1] + '.stride')):

                with open(outPath + '/dssp/unbound/' + targetName + '_' + complexDir[c].split('_')[1] + '/' +
                              files_rr[f].split('.rr')[0] + '_' + complexDir[c].split('_')[1] + '.stride') as stride1File:
                    for line in stride1File:
                        tmp = line.split()
                        if(len(tmp) > 0 and line[0:(0+3)] == "ASG"):
                            chain1Stride.append(line)

            if(os.path.exists(outPath + '/dssp/unbound/' + targetName + '_' + complexDir[c].split('_')[2] + '/' +
                               files_rr[f].split('.rr')[0] + '_' + complexDir[c].split('_')[2] + '.dssp')):

                with open(outPath + '/dssp/unbound/' + targetName + '_' + complexDir[c].split('_')[2] + '/' +
                              files_rr[f].split('.rr')[0] + '_' + complexDir[c].split('_')[2] + '.dssp') as dssp2File:
                    for line in dssp2File:
                        if (line2Count<1):
                            if (line[2:(2+1)] == '#'):
                                line2Count += 1
                                continue
                        if(line2Count > 0):
                            if(len(line) > 0):
                                chain2Dssp.append(line)
                    dsspChain2Found = 1
                    
            elif(os.path.exists(outPath + '/dssp/unbound/' + targetName + '_' + complexDir[c].split('_')[2] + '/' +
                              files_rr[f].split('.rr')[0] + '_' + complexDir[c].split('_')[2] + '.stride')):
                with open(outPath + '/dssp/unbound/' + targetName + '_' + complexDir[c].split('_')[2] + '/' +
                              files_rr[f].split('.rr')[0] + '_' + complexDir[c].split('_')[2] + '.stride') as stride2File:
                    for line in stride2File:
                        tmp = line.split()
                        if(len(tmp) > 0 and line[0:(0+3)] == "ASG"):
                            chain2Stride.append(line)
                            
            #------get sequence length-------#
            #                                #
            #--------------------------------#
            chain1SeqLen = fastaLength(outPath + '/fasta/unbound/' + targetName + '_' + complexDir[c].split('_')[1] + '.fasta')
            chain2SeqLen = fastaLength(outPath + '/fasta/unbound/' + targetName + '_' + complexDir[c].split('_')[2] + '.fasta')

            neffChain1 = ''
            neffChain2 = ''
            neffChainAll = ''

            if(os.path.exists(outPath + '/neff/unbound/' + targetName + '_' + complexDir[c].split('_')[1] + '.neff')):
                neffChain1 = get_neff(outPath + '/neff/unbound/' + targetName + '_' + complexDir[c].split('_')[1] + '.neff')
                
            if(os.path.exists(outPath + '/neff/unbound/' + targetName + '_' + complexDir[c].split('_')[2] + '.neff')):
                neffChain2 = get_neff(outPath + '/neff/unbound/' + targetName + '_' + complexDir[c].split('_')[2] + '.neff')

            if(os.path.exists(outPath + '/neff/bound/' + targetName + '_' + complexDir[c].split('_')[1] + '_' + complexDir[c].split('_')[2] + '.neff')):
                neffChainAll = get_neff(outPath + '/neff/bound/' + targetName + '_' + complexDir[c].split('_')[1] + '_' + complexDir[c].split('_')[2] + '.neff')

            #------get angle features------#
            #                              #
            #------------------------------#
            if(os.path.exists(outPath + '/orientation/' + targetName + '_' + complexDir[c].split('_')[1] + '_' + complexDir[c].split('_')[2] + '/' +
                              files_rr[f].split('.rr')[0] + '.ori')):
                with open(outPath + '/orientation/' + targetName + '_' + complexDir[c].split('_')[1] + '_' + complexDir[c].split('_')[2] + '/' +
                              files_rr[f].split('.rr')[0] + '.ori') as oFile:
                    orientLines = []
                    for line in oFile:
                       orientLines.append(line)


            #---generate feat for pairs-----#
            #                               #
            #-------------------------------#
            outFile_feat = open(outPath + '/features/' + complexDir[c] + '/' + targetName + '_' + files_rr[f].split('.rr')[0] + '.feat', 'w')
            for d in range(len(decLines)):
                try:
                    tmpD = decLines[d].split()
                    aaFeat1 = ''
                    aaFeat2 = ''
                    ss1 = ''
                    sa1 = ''
                    aa1 = ''
                    aa2 = ''
                    
                    ss2 = ''
                    sa2 = ''
                    tco1 = ''
                    tco2 = ''
                    kappa1 = ''
                    kappa2 = ''
                    alpha1 = ''
                    alpha2 = ''
                    phi1 = ''
                    phi2 = ''
                    psi1 = ''
                    psi2 = ''
                    sinPhi1 = ''
                    cosPhi1 = ''
                    sinPsi2 = ''
                    cosPsi2 = ''
                    
                    ssFeat1 = ''
                    ssFeat2 = ''
                    saFeat1 = ''
                    saFeat2 = ''
                    tcoFeat = ''
                    kappaFeatSin = ''
                    alphaFeatSin = ''
                    phiFeatSin = ''
                    psiFeatSin = ''

                    kappaFeatCos = ''
                    alphaFeatCos = ''
                    phiFeatCos = ''
                    psiFeatCos = ''

                    sinPhi = ''
                    sinOmega = ''
                    sinTheta = ''
                    cosPhi = ''
                    cosOmega = ''
                    cosTheta = ''

                    energFeat = ''

                    pssmUnbound1 = ''
                    pssmUnbound2 = ''

                    pssmBound1 = ''
                    pssmBound2 = ''

                    feat = ''

                    #dssp
                    if(len(chain1Dssp) > 0):
                        for ds1 in range(len(chain1Dssp)):
                            if(tmpD[0] == chain1Dssp[ds1][6:(6+4)].strip()):


                                relResFeat1 = int(chain1Dssp[ds1].split()[0]) / chain1SeqLen
                                aa1 = chain1Dssp[ds1][13:(13+1)]
                        
                                ss1 = get8to3ss(chain1Dssp[ds1][16:(16+1)])
                                if(isfloat(chain1Dssp[ds1][35:(35+3)].strip())):
                                    sa1 = float(chain1Dssp[ds1][35:(35+3)])

                                if(isfloat(chain1Dssp[ds1][103:(103+6)].strip())):
                                    phi1 = float(chain1Dssp[ds1][103:(103+6)].strip())
                                if(isfloat(chain1Dssp[ds1][109:(109+6)].strip())):
                                    psi1 = float(chain1Dssp[ds1][109:(109+6)].strip())

                                sinPhi1 = sigmoid(math.sin(phi1))
                                cosPhi1 = sigmoid(math.cos(phi1))

                                sinPsi1 = sigmoid(math.sin(psi1))
                                cosPsi1 = sigmoid(math.cos(psi1))
                            
                                break
                            
                    elif(len(chain1Stride) > 0):
                        for std1 in range(len(chain1Stride)):
                            if(tmpD[0] == chain1Stride[std1][11:(11+4)].strip()):


                                relResFeat1 = int(chain1Stride[std1][16:(16+4)].strip()) / chain1SeqLen #extract the serial number, not the real actual sequence number in the pdb

                                aa1 = get3to1aa(chain1Stride[std1][5:(5+3)].strip())
                                
                                ss1 = chain1Stride[std1][24:(24+1)]
                                if(isfloat(chain1Stride[std1][61:(61+8)].strip())):
                                    sa1 = float(chain1Stride[std1][61:(61+8)].strip())

                                if(isfloat(chain1Stride[std1][42:(42+7)].strip())):
                                    phi1 = float(chain1Stride[std1][42:(42+7)].strip())
                                if(isfloat(chain1Stride[std1][52:(52+7)].strip())):
                                    psi1 = float(chain1Stride[std1][52:(52+7)].strip())

                                sinPhi1 = sigmoid(math.sin(phi1))
                                cosPhi1 = sigmoid(math.cos(phi1))

                                sinPsi1 = sigmoid(math.sin(psi1))
                                cosPsi1 = sigmoid(math.cos(psi1))
                                break
                                


            
                    if(len(chain2Dssp) > 0):
                        for ds2 in range(len(chain2Dssp)):
                            if(tmpD[1] == chain2Dssp[ds2][6:(6+4)].strip()):

                                relResFeat2 = int(chain2Dssp[ds2].split()[0]) / chain2SeqLen
                                
                                aa2 = chain2Dssp[ds2][13:(13+1)]
                                
                                ss2 = get8to3ss(chain2Dssp[ds2][16:(16+1)])
                                if(isfloat(chain2Dssp[ds2][35:(35+3)])):
                                    sa2 = float(chain2Dssp[ds2][35:(35+3)])

                                if(isfloat(chain2Dssp[ds2][103:(103+6)])):
                                    phi2 = float(chain2Dssp[ds2][103:(103+6)])
                                if(isfloat(chain2Dssp[ds2][109:(109+6)])):
                                    psi2 = float(chain2Dssp[ds2][109:(109+6)])

                                sinPhi2 = sigmoid(math.sin(phi2))
                                cosPhi2 = sigmoid(math.cos(phi2))

                                sinPsi2 = sigmoid(math.sin(psi2))
                                cosPsi2 = sigmoid(math.cos(psi2))

                                break
                            
                    elif(len(chain2Stride) > 0):
                        for std2 in range(len(chain2Stride)):
                            if(tmpD[1] == chain2Stride[std2][11:(11+4)].strip()):


                                relResFeat2 = int(chain2Stride[std2][16:(16+4)].strip()) / chain2SeqLen #extract the serial number, not the real actual sequence number in the pdb

                                aa2 = get3to1aa(chain2Stride[std2][5:(5+3)].strip())
                                
                                ss2 = chain2Stride[std2][24:(24+1)]
                                if(isfloat(chain2Stride[std2][61:(61+8)].strip())):
                                    sa2 = float(chain2Stride[std2][61:(61+8)].strip())

                                if(isfloat(chain2Stride[std2][42:(42+7)].strip())):
                                    phi2 = float(chain2Stride[std2][42:(42+7)].strip())
                                if(isfloat(chain2Stride[std2][52:(52+7)].strip())):
                                    psi2 = float(chain2Stride[std2][52:(52+7)].strip())

                                sinPhi2 = sigmoid(math.sin(phi2))
                                cosPhi2 = sigmoid(math.cos(phi2))

                                sinPsi2 = sigmoid(math.sin(psi2))
                                cosPsi2 = sigmoid(math.cos(psi2))
                                
                                break

                    aaFeat1 = aaGroup(aa1)
                    aaFeat2 = aaGroup(aa2)
                    
                    ssFeat1 = get8to3ssOneHot(ss1)
                    ssFeat2 = get8to3ssOneHot(ss2)


                    #orientation
                    for o in range(len(orientLines)):
                        tmpO = orientLines[o].split()
                        if(int(tmpO[0]) == int(tmpD[0]) and int(tmpO[1]) == int(tmpD[1])):
                            sinPhi1_or = sigmoid(math.sin(float(tmpO[2])))
                            cosPhi1_or = sigmoid(math.cos(float(tmpO[2])))
                            sinPhi2_or = sigmoid(math.sin(float(tmpO[3])))
                            cosPhi2_or = sigmoid(math.cos(float(tmpO[3])))
                            sinTheta1_or = sigmoid(math.sin(float(tmpO[4])))
                            cosTheta1_or = sigmoid(math.cos(float(tmpO[4])))
                            sinTheta2_or = sigmoid(math.sin(float(tmpO[5])))
                            cosTheta2_or = sigmoid(math.cos(float(tmpO[5])))
                            sinOmega_or = sigmoid(math.sin(float(tmpO[6])))
                            cosOmega_or = sigmoid(math.cos(float(tmpO[6])))

                            break

                    saFeat1 = getSolAccOneHot(aa1, float(sa1))
                    saFeat2 = getSolAccOneHot(aa2, float(sa2))

                    phiFeatSin = sigmoid((math.sin(phi1) - math.sin(phi2)) **2)
                    psiFeatSin = sigmoid((math.sin(psi1) - math.sin(psi2)) **2)

                    phiFeatCos = sigmoid((math.cos(phi1) - math.cos(phi2)) **2)
                    psiFeatCos = sigmoid((math.cos(psi1) - math.cos(psi2)) **2)

                    featLeft = str(neffChain1) + ' ' + str(neffChainAll) + ' ' + aaFeat1 + ' ' + ssFeat1 + ' ' + \
                               saFeat1 + ' ' + str(relResFeat1) + ' ' + str(sinPhi1) + ' ' + str(cosPhi1) + ' ' + \
                               str(sinPsi1) + ' ' + str(cosPsi1)
                    
                    featRight= str(neffChain2) + ' ' + str(neffChainAll) + ' ' + aaFeat2 + ' ' + ssFeat2 + ' ' + saFeat2 + ' ' + \
                               str(relResFeat2) + ' ' + str(sinPhi2) + ' ' + str(cosPhi2) + ' ' + str(sinPsi2) + ' ' + str(cosPsi2)
                    
                    edgeFeat = str(sinPhi1_or) + ' ' + str(cosPhi1_or) + ' ' + str(sinPhi2_or) + ' ' + str(cosPhi2_or) + ' ' + str(sinTheta1_or) + ' ' + \
                               str(cosTheta1_or) + ' ' + str(sinTheta2_or) + ' ' + str(cosTheta2_or) + ' ' + str(sinOmega_or) + ' ' + str(cosOmega_or) + ' ' + \
                               edgeFeatEncoding(float(tmpD[2]))

                    if(len(featLeft.split()) == 17 and len(featRight.split()) == 17 and len(edgeFeat.split()) == 27):
                        feat = tmpD[0] + ' ' + tmpD[1] + ',' + featLeft + ',' + featRight + ',' + edgeFeat
                        outFile_feat.write(feat + '\n')
                
                except Exception as e:
                    print("Error occurred for: " + targetName + ' ' + files_rr[f] + ' ' + str(e))


            outFile_feat.close()

def calNumInterface(complex_name, pdb):
    total = 0
    with open(outPath + '/distance/' + complex_name + '/' + pdb.split('.pdb')[0] + '.rr') as dFile:
        for line in dFile:
            if(float(line.split()[2]) < intFDist):
                total += 1
    return total

def renumber(iterable):
    seen = {}
    counter = 0
    renum_list = []
    ori_list = []

    for x in iterable:
        i = seen.get(x)

        if i is None:
            seen[x] = counter
            renum_list.append(counter)
            counter += 1
        else:
            renum_list.append(i)
        ori_list.append(x)
    return renum_list, ori_list

class GATLayer(nn.Module):
    def __init__(self, in_dim, out_dim, efeats):
        super(GATLayer, self).__init__()

        self.W_msg = nn.Linear(in_dim + efeats, out_dim)
        self.W_apply = nn.Linear(in_dim + out_dim, out_dim)
        self.reset_parameters()

    def reset_parameters(self):
        gain = nn.init.calculate_gain('relu')
        nn.init.xavier_normal_(self.W_msg.weight, gain=gain)
        nn.init.xavier_normal_(self.W_apply.weight, gain=gain)

    def aggr(self, msg):
        return fn.sum(msg, 'msg_attn')

    def msg(self, int_edges):

        msg_edge = torch.cat([int_edges.src['h'], int_edges.data['h']], dim = 1)
        msg_edge_f = self.W_msg(msg_edge)
        return {'msg': F.relu(msg_edge_f)}

    def forward(self, g_dgl, node_feats, edge_feats):
        with g_dgl.local_scope():
            g = g_dgl
            g.ndata['h'] = node_feats
            g.edata['h'] = edge_feats
            g.update_all(self.msg, self.aggr('msg'))
            h_dat = [g.ndata['h'], g.ndata['msg_attn']]
            z = torch.cat(h_dat, dim = 1)
            a = F.relu(self.W_apply(z))
            g.ndata['h_'] = a
            #g.update_all(fn.copy_u('h_', 'r'), fn.sum('r', 'h'))
            return g.ndata['h_']

class MultiHeadGATLayer(nn.Module):
    def __init__(self, in_dim, out_dim, num_heads, efeats, merge='cat'):
        super(MultiHeadGATLayer, self).__init__()
        self.heads = nn.ModuleList()
        for i in range(num_heads):
            self.heads.append(GATLayer(in_dim, out_dim, efeats))
        self.merge = merge

    def forward(self, g, h, efeats):
        head_outs = [attn_head(g, h, efeats) for attn_head in self.heads]
        if self.merge == 'cat':
            return torch.cat(head_outs, dim=1)
        else:
            return torch.mean(torch.stack(head_outs))

class GAT(nn.Module):
    def __init__(self, in_dim, hidden_dim, out_dim, num_heads, efeats, num_layers):
        super(GAT, self).__init__()

        self.layers = nn.ModuleList()

        #input layer
        self.layers.append(MultiHeadGATLayer(in_dim, hidden_dim, num_heads, efeats))
        #add layer recursively
        for i in range(num_layers):
            hidden_dim_in = hidden_dim
            hidden_dim_out = int(hidden_dim_in * 0.5)
            self.layers.append(MultiHeadGATLayer(hidden_dim_in * num_heads, hidden_dim_out, num_heads, efeats))
            hidden_dim = hidden_dim_out
        #add output layer
        self.layers.append(MultiHeadGATLayer(hidden_dim_out * num_heads, out_dim, 1, efeats))
        self.dropout = nn.Dropout(p=0.5)
        
    def forward(self, g, h, efeats):
        for i, layer in enumerate(self.layers):
            h = self.dropout(h)            
            h = layer(g, h, efeats)
        g.ndata['h'] = h
        g.apply_edges(fn.u_dot_v('h', 'h', 'score'))
        return g.edata['score'].squeeze(1)

def score():
    print("Scoring...")
    os.chdir(working_path) #just in case
    featurePath = os.listdir(outPath + '/features/')
    #for each complex
    for f in range(len(featurePath)):
        if not(os.path.exists(outPath + '/scores/' + featurePath[f])):
            os.system('mkdir -p ' + outPath + '/scores/' + featurePath[f])

        saved_model = PIQLE_path + 'model/model_weight'

        #-----------------make prediction--------------#
        #                                              #
        #----------------------------------------------#
        featFiles = os.listdir(outPath + '/features/' + featurePath[f])

        for j in range(len(featFiles)):
            featureFiles = featFiles[j]
            distLines = []

            nodes1 = []
            nodes2 = []

            nodes1FeatTrack = {}
            nodes2FeatTrack = {}
            nodes1Feat = []
            nodes2Feat = []
            labels1 = []
            labels2 = []
            
            nodeFeats = []
            edgeFeats = []
            labels_all = []
            c1 = 0
            c2 = 0
            if((outPath + '/features/' + featurePath[f] + '/' + featureFiles).endswith('.feat') and
               os.stat(outPath + '/features/' + featurePath[f] + '/' + featureFiles).st_size != 0):
                print(outPath + '/features/' + featurePath[f] + '/' + featureFiles)
                with open(outPath + '/features/' + featurePath[f] + '/' + featureFiles) as fFile:
                    for line in fFile:
                        tmp = line.split(',')

                        residuePair = tmp[0].split()
                        leftFeatInfo = tmp[1].split()
                        rightFeatInfo = tmp[2].split()
                        edgeFeatInfo = tmp[3].split()
                        #labelInfo = tmp[4] #only one value, so no split function applied
                        
                        nodes1.append(int(residuePair[0]))
                        nodes2.append(int(residuePair[1]))

                        if(nodes1FeatTrack.get(int(residuePair[0])) == None):
                            nodes1FeatTrack[int(residuePair[0])] = c1
                            nodes1Feat.append([float(i) for i in leftFeatInfo])
                            c1 += 1
                        if(nodes2FeatTrack.get(int(residuePair[1])) == None):
                            nodes2FeatTrack[int(residuePair[1])] = c2
                            nodes2Feat.append([float(i) for i in rightFeatInfo])
                            c2 += 1

                        edgeFeats.append([float(i) for i in edgeFeatInfo])
                        #labels_all.append(float(labelInfo))


                #process node features to append sequentially
                for n1 in range(len(nodes1Feat)):
                    nodeFeats.append(nodes1Feat[n1])
                for n2 in range(len(nodes2Feat)):
                    nodeFeats.append(nodes2Feat[n2])
                
                #renumber the nodes starting from 0
                nodes1 = renumber(nodes1)[0]
                #process nodes 2:
                #for the second list, add last element of the first list for continuation
                #(see the paper prototype)
                nodes2 = renumber(nodes2)[0]
                nodes2 = list([x + nodes1[-1] + 1 for x in nodes2]) #add the last element of list 1 to every element of list 2

                #print(len(nodes1Feat))
                #print(len(nodes2Feat))
            
                #define the graph
                nodes1 = th.tensor(nodes1)
                nodes2 = th.tensor(nodes2)
                graph = dgl.graph((nodes1, nodes2))

                #process the features and label
                nodeFeats =  th.tensor(nodeFeats)
                edgeFeats =  th.tensor(edgeFeats)
                #labels_all = th.tensor(labels_all)

                #add features and label to the graph
                graph.ndata['n'] = nodeFeats
                graph.edata['e'] = edgeFeats
                #graph.edata['l'] = labels_all

                #load model
                model = th.load(saved_model, map_location=torch.device('cpu'))
                model.eval()
                nfeats = graph.ndata['n']
                efeats = graph.edata['e']
                output = model(graph, nfeats, efeats)
                #print("Prediction output")
                outFile = open(outPath + '/scores/' + featurePath[f] + '/' + featureFiles.split('.feat')[0] + '_pred', 'w')
                #labels_all = labels_all.tolist()
                #print(labels_all[0])
                for i in range(len(output)):
                    #outFile.write(featureFiles.split('.')[0] + ' ' + str(output[i].detach().numpy()) + ' ' + str(labels_all[i]) + '\n')
                    outFile.write(featureFiles.split('.')[0] + ' ' + str(output[i].detach().numpy()) + '\n')

def finalizeScore():
    outFile = open(outPath + '/' + targetName + '.PIQLE', 'w')

    decoys = os.listdir(decoyPath)
    scoreDirs = os.listdir(outPath + '/scores')

    #header line
    header_line = ''
    for n in range(len(scoreDirs)):
        header_line += scoreDirs[n] + ' '

    #calculate global score for each decoy
    for d in range(len(decoys)):
        #look for all the score directories
        global_score_weighted = 0
        global_score_average = 0
        scores_line = ''
        indiv_score_weighted = []
        indiv_score = []
        total_weight = 0
        
        for s in range(len(scoreDirs)):
            scoreFiles = os.listdir(outPath + '/scores/' + scoreDirs[s])
            for f in range(len(scoreFiles)):
                local_scores = []
                if(targetName + '_' + decoys[d].split('.pdb')[0] == scoreFiles[f].split('_pred')[0]):
                    with open(outPath + '/scores/' + scoreDirs[s] + '/' + scoreFiles[f]) as sFile:
                        for line in sFile:
                            tmp = line.split()
                            local_scores.append(float(tmp[1]))

                    weight = calNumInterface(scoreDirs[s], decoys[d])
                    total_weight += weight
                    scores_line += str(sum(local_scores) / (len(local_scores))) + ' ' + str(weight) + ' '
                    indiv_score.append(sum(local_scores) / (len(local_scores)))
                    indiv_score_weighted.append(sum(local_scores) / (len(local_scores)) * weight)
                    break

        if(len(indiv_score) > 0):
            global_score_average = sum(indiv_score) / len(indiv_score)
            outFile.write(decoys[d] + ' ' + str(global_score_average) + '\n')
    
    outFile.close()
    #sort the score file#
    with open(outPath + '/' + targetName + '.PIQLE') as sFile:
        lines = []
        for line in sFile:
                tmp = line.split()
                if(len(tmp) > 0):
                        lines.append(line)

    score_file = open(outPath + '/' + targetName + '.PIQLE', 'w')
    for line in sorted(lines, key=lambda line: float(line.split()[1]), reverse = True):
        score_file.write(line)
    score_file.close()

    if(os.path.exists(outPath + '/' + targetName + '.PIQLE') and os.path.getsize(outPath + '/' + targetName + '.PIQLE') > 0):
            print("Congratulations! All processes are done successfully")

def main():      
    runDSSP()
    generatePairs()
    generateOrientation()
    generateDistance()
    concatMSA()
    calculateNeff()
    generateFeatures()
    score()
    finalizeScore()
        
if __name__ == '__main__':
        main()
