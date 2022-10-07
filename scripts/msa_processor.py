import os, sys, fnmatch, math
import numpy as np
import os,sys,shutil, re, subprocess, optparse, multiprocessing

parser = optparse.OptionParser()

parser.add_option('--aln', dest='alnDir',
        default = '',    # default empty!
        help = 'Alignment directory')
parser.add_option('--out', dest='outPath',
        default = '',    # default empty!
        help = 'Output path')

(options,args) = parser.parse_args()

alnLoc = options.alnDir
outPath = options.outPath

files = os.listdir(alnLoc)

for i in range(len(files)):
	if(files[i].endswith('.concat.aln')):
		lines = []
		with open(alnLoc + '/' + files[i]) as aFile:
			for line in aFile:
				if(line[0] != '>'):
					lines.append(line)
		outFile = open(outPath + '/' + files[i], 'w')
		for j in range(len(lines)):
			outFile.write(lines[j])
