#from scipy import optimize
#from scipy.optimize import minimize
import os,math
import optparse    # for option sorting

parser = optparse.OptionParser()
parser.add_option('-d', dest='decoyDir',
        default = '',    # default empty!
        help = 'Decoy directory')

parser.add_option('-l', dest='chain1',
        default = '',    # default empty!
        help = 'First chain')

parser.add_option('-r', dest='chain2',
        default = '',    # default empty!
        help = 'Second chain')

parser.add_option('-o', dest='outputDir',
        default = '',    # default empty!
        help = 'Output directory')

(options,args) = parser.parse_args()

decoyDir = options.decoyDir
chain1 = options.chain1
chain2 = options.chain2
outPut = options.outputDir

def get3to1aa(aa):
    dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    return dict[aa]

pdbList = os.listdir(decoyDir)
for d in range(len(pdbList)):
    print(decoyDir + '/' + pdbList[d])
    
    pdbFile = decoyDir + '/' + pdbList[d]

    chainRecord = []
    lines = []

    with open(pdbFile) as pFile:
        for line in pFile:
            if(len(line.split())>0 and line.split()[0] == "ATOM"):
                chainRecord.append(line[21:(21+1)].strip())
                lines.append(line.strip())
                
    uniqueChain = list(set(chainRecord))
    chain1Res = []
    chain1ResNum = []
    chain1CbLine = []
    chain2Res = []
    chain2ResNum = []
    chain2CbLine= []
    
    if(len(uniqueChain) == 2): #process only if 2 chains
        #chain1
        for m in range(len(lines)):
            if(lines[m][12:(12+4)].strip() == 'CB' or (lines[m][17:(17+3)] == "GLY" and lines[m][12:(12+4)].strip() == 'CA')):
                if(lines[m][21:(21+1)] == chain1):
                    #chain1Res.append(get3to1aa(lines[m][17:(17+3)]))
                    chain1CbLine.append(lines[m])
                    chain1ResNum.append(lines[m][22:(22+4)])
        #chain2
        for m in range(len(lines)):
            if(lines[m][12:(12+4)].strip() == 'CB' or (lines[m][17:(17+3)] == "GLY" and lines[m][12:(12+4)].strip() == 'CA')):
                if(lines[m][21:(21+1)] == chain2):
                    #chain2Res.append(get3to1aa(lines[m][17:(17+3)]))
                    chain2CbLine.append(lines[m])
                    chain2ResNum.append(lines[m][22:(22+4)])
        outFile = open(outPut + '/' + pdbList[d].split('.pdb')[0] + '.rr', 'w')
        for c1 in range(len(chain1ResNum)):
            for c2 in range(len(chain2ResNum)):
                Xi = float(chain1CbLine[c1][30:(30+8)].strip())
                Yi = float(chain1CbLine[c1][38:(38+8)].strip())
                Zi = float(chain1CbLine[c1][46:(46+8)].strip())

                Xj = float(chain2CbLine[c2][30:(30+8)].strip())
                Yj = float(chain2CbLine[c2][38:(38+8)].strip())
                Zj = float(chain2CbLine[c2][46:(46+8)].strip())
                           
                dist = math.sqrt((Xi - Xj) ** 2 + (Yi - Yj) ** 2 + (Zi - Zj) ** 2)
                #print(str(chain1ResNum[c1].strip()) + ' ' + str(chain2ResNum[c2].strip()) + ' ' + str(dist) + '\n')
                outFile.write(str(chain1ResNum[c1].strip()) + ' ' + str(chain2ResNum[c2].strip()) + ' ' + str(dist) + ' ' + 
                              chain1CbLine[c1][21:(21+1)].strip() + ' ' + chain2CbLine[c2][21:(21+1)].strip() + '\n')

            
    
    
            
