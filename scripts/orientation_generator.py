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

def get_distance(x1, x2, x3, y1, y2, y3):
        return math.sqrt((x1 - y1) ** 2 + (x2 - y2) ** 2 + (x3 - y3) ** 2)

def getAngle(x1, x2, x3, y1, y2, y3, z1, z2, z3):
        acc = 0.0
        d1 = 0
        d2 = 0
        acc = (x2 - x1) * (x2 - x3) + (y2 - y1) * (y2 - y3) + (z2 - z1) * (z2 - z3)
        d1 = math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
        d2 = math.sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
        if (d1 == 0 or d2 == 0):
                return 0
        acc = acc / (d1 * d2);
        if (acc > 1.0):
                acc = 1.0
        elif (acc < -1.0):
                acc = -1.0;
        acc = math.acos(acc);
        return acc;


def getDihedral(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4):

        #getDifferences
        #d1 = difference(p1, p2)
        qx = x2 - x1;
        qy = y2 - y1;
        qz = z2 - z1;

        #d2 = difference(p3, p2)
        rx = x2 - x3;
        ry = y2 - y3;
        rz = z2 - z3;

        #d3 = difference(p4, p3)
        sx = x3 - x4;
        sy = y3 - y4;
        sz = z3 - z4;

        #get cross product#
        #cp1 = crossProduct(d1, d2)
        tx = qy * rz - qz * ry;
        ty = qz * rx - qx * rz;
        tz = qx * ry - qy * rx;

        #cp2 = crossProduct(d3, d1)

        ux = sy * rz - sz * ry;
        uy = sz * rx - sx * rz;
        uz = sx * ry - sy * rx;

        #cp3 = crossProduct(cp2, cp1)
        vx = uy * tz - uz * ty;
        vy = uz * tx - ux * tz;
        vz = ux * ty - uy * tx;

        #get dot product#
        #dot = getDotProduct(cp3, d2)
        w = vx * rx + vy * ry + vz * rz;

        #get angle#
        zx = zy = zz = 0.0;
        accTau = 0.0;
        accTau = (zx - tx) * (zx - ux) + (zy - ty) * (zy - uy) + (zz - tz) * (zz - uz);
        d1Tau = math.sqrt((tx - zx) * (tx - zx) + (ty - zy) * (ty - zy) + (tz - zz) * (tz - zz));
        d2Tau = math.sqrt((ux - zx) * (ux - zx) + (uy - zy) * (uy - zy) + (uz - zz) * (uz - zz));

        if (d1Tau == 0 or d2Tau == 0):
                return 0

        accTau = accTau / (d1Tau * d2Tau);

        if (accTau > 1.0):
                accTau = 1.0;
        elif (accTau < -1.0):
                accTau = -1.0;
        accTau = math.acos(accTau);
        #end of getAngle#
        if (w < 0):
                accTau = -accTau;
        return accTau

#showing angle (for future use)                                                                                                                                                                                                                                  
def showAngles(pos):
        x = 0
        y = 1
        z = 2
        #f = open('penaltya.cst', 'w')
        for i in range(1, len(pos)-2):
                print ('res no ' + str(i) + ', secondary: ' + str(secondary[0][i]))
                print ('dihedral angle : ' + str(180 / 3.1416 * getDihedral(pos[i-1][x], pos[i][x], pos[i+1][x], pos[i+2][x], pos[i-1][y], pos[i][y], pos[i+1][y], pos[i+2][y], pos[i-1][z], pos[i][z], pos[i+1][z], pos[i+2][z])))
                print ('tau angle : ' +str(180 / 3.1416 * getAngle(pos[i-1][x], pos[i][x], pos[i+1][x], pos[i-1][y], pos[i][y], pos[i+1][y], pos[i-1][z], pos[i][z], pos[i+1][z])))
                #f.write('Angle CA ' + str(position) + ' CA ' + str(position) + ' LINEAR_PENALTY ' + str(d) + ' 0 0 1.0\n')


def cal_orientation(pdb, outFile, chainIds):

        #---check the chains---#
        #                      #
        #----------------------#
        chainRecord = []
        with open(pdb) as pFile:
            for line in pFile:
                if(len(line.split())>0 and line.split()[0] == "ATOM"):
                    chainRecord.append(line[21:(21+1)].strip())
        uniqueChain = list(set(chainRecord))

        if(len(uniqueChain) == 2): #process only if dimer
                fpdb = open(pdb, 'r')
                #pdb_atoms = fpdb.readlines()
                flines = fpdb.readlines()
                pos1 = []
                pos2 = []                
                for atomline in flines: 
                        if(atomline[:4] == 'ATOM'):  
                                if(atomline[21] == chainIds[0]):
                                        pos1.append(atomline) #[atomline[21], atomline[22:26].strip(), float(atomline[30:38].strip()), float(atomline[38:46].strip()), float(atomline[46:54].strip())])     
                                elif(atomline[21] == chainIds[1]):  
                                        pos2.append(atomline) #[21], atomline[22:26].strip(), float(atomline[30:38].strip()), float(atomline[38:46].strip()), float(atomline[46:54].strip())]) 
                fpdb.close()   
                             
                Ca_info = {}
                Cb_info = {}
                N_info = {}
                for line in pos1:
                        if(line[:4] == "ATOM" and line[17:20].strip() != "GLY"):
                                  
                              
                                Ca = []
                                Cb = []
                                N = []
                                res_no = 0
                                if(line[12:16].strip() == "CA"):
                                        x = float(line[30:38].strip())
                                        y = float(line[38:46].strip())
                                        z = float(line[46:54].strip())
                                        Ca = [x, y, z]
                                        res_no = int(line[22:26].strip())
                                        Ca_info[res_no] = Ca

                                if(line[12:16].strip() == "CB"):
                                        x = float(line[30:38].strip())
                                        y = float(line[38:46].strip())
                                        z = float(line[46:54].strip())
                                        Cb = [x, y, z]
                                        res_no = int(line[22:26].strip())
                                        Cb_info[res_no] = Cb

                                if(line[12:16].strip() == "N"):
                                        x = float(line[30:38].strip())
                                        y = float(line[38:46].strip())
                                        z = float(line[46:54].strip())
                                        N = [x, y, z]
                                        res_no = int(line[22:26].strip())
                                        N_info[res_no] = N

                        if(line[:4] == "ATOM" and line[17:20].strip() == "GLY"):
                                  
                              
                                Ca = []
                                Cb = []
                                N = []
                                res_no = 0
                                if(line[12:16].strip() == "CA"):
                                        x = float(line[30:38].strip())
                                        y = float(line[38:46].strip())
                                        z = float(line[46:54].strip())
                                        Ca = [x, y, z]
                                        res_no = int(line[22:26].strip())
                                        Ca_info[res_no] = Ca

                                        Cb_info[res_no] = Ca

                                if(line[12:16].strip() == "N"):
                                        x = float(line[30:38].strip())
                                        y = float(line[38:46].strip())
                                        z = float(line[46:54].strip())
                                        N = [x, y, z]
                                        res_no = int(line[22:26].strip())
                                        N_info[res_no] = N


                          

                Ca_info2 = {}
                Cb_info2 = {}
                N_info2 = {}
                for line in pos2:
                        if(line[:4] == "ATOM" and line[17:20].strip() != "GLY"):


                                Ca = []
                                Cb = []
                                N = []
                                res_no = 0
                                if(line[12:16].strip() == "CA"):
                                        x = float(line[30:38].strip())
                                        y = float(line[38:46].strip())
                                        z = float(line[46:54].strip())
                                        Ca = [x, y, z]
                                        res_no = int(line[22:26].strip())
                                        Ca_info2[res_no] = Ca

                                if(line[12:16].strip() == "CB"):
                                        x = float(line[30:38].strip())
                                        y = float(line[38:46].strip())
                                        z = float(line[46:54].strip())
                                        Cb = [x, y, z]
                                        res_no = int(line[22:26].strip())
                                        Cb_info2[res_no] = Cb

                                if(line[12:16].strip() == "N"):
                                        x = float(line[30:38].strip())
                                        y = float(line[38:46].strip())
                                        z = float(line[46:54].strip())
                                        N = [x, y, z]
                                        res_no = int(line[22:26].strip())
                                        N_info2[res_no] = N

                        if(line[:4] == "ATOM" and line[17:20].strip() == "GLY"):


                                Ca = []
                                Cb = []
                                N = []
                                res_no = 0
                                if(line[12:16].strip() == "CA"):
                                        x = float(line[30:38].strip())
                                        y = float(line[38:46].strip())
                                        z = float(line[46:54].strip())
                                        Ca = [x, y, z]
                                        res_no = int(line[22:26].strip())
                                        Ca_info2[res_no] = Ca

                                        Cb_info2[res_no] = Ca

                                if(line[12:16].strip() == "N"):
                                        x = float(line[30:38].strip())
                                        y = float(line[38:46].strip())
                                        z = float(line[46:54].strip())
                                        N = [x, y, z]
                                        res_no = int(line[22:26].strip())
                                        N_info2[res_no] = N

                print(len(Ca_info2))
                print(len(Cb_info2))
                for res_no in Ca_info:
                        for res_no2 in Ca_info2:
                                try:
                                        cb_cb_distance = get_distance(Cb_info[res_no][0], Cb_info[res_no][1],  Cb_info[res_no][2], Cb_info2[res_no2][0], Cb_info2[res_no2][1], Cb_info2[res_no2][2])

                                        #if(res_no != res_no2): #does not apply for dimer
                                        phi = getAngle(Ca_info[res_no][0], Cb_info[res_no][0], Cb_info2[res_no2][0],
                                                       Ca_info[res_no][1], Cb_info[res_no][1], Cb_info2[res_no2][1],
                                                       Ca_info[res_no][2], Cb_info[res_no][2], Cb_info2[res_no2][2])

                                        phi2 = getAngle(Ca_info2[res_no2][0], Cb_info2[res_no2][0], Cb_info[res_no][0],
                                                       Ca_info2[res_no2][1], Cb_info2[res_no2][1], Cb_info[res_no][1],
                                                       Ca_info2[res_no2][2], Cb_info2[res_no2][2], Cb_info[res_no][2])
                                        
                                        #print(str(res_no) + "," + str(res_no2) + " omega = " + str(omega))
                                        theta = getDihedral(N_info[res_no][0], Ca_info[res_no][0], Cb_info[res_no][0], Cb_info2[res_no2][0],
                                                            N_info[res_no][1], Ca_info[res_no][1], Cb_info[res_no][1], Cb_info2[res_no2][1],
                                                            N_info[res_no][2], Ca_info[res_no][2], Cb_info[res_no][2], Cb_info2[res_no2][2])

                                        #print(str(res_no) + "," + str(res_no2) + " omega = " + str(omega))
                                        theta2 = getDihedral(N_info2[res_no2][0], Ca_info2[res_no2][0], Cb_info2[res_no2][0], Cb_info[res_no][0],
                                                             N_info2[res_no2][1], Ca_info2[res_no2][1], Cb_info2[res_no2][1], Cb_info[res_no][1],
                                                             N_info2[res_no2][2], Ca_info2[res_no2][2], Cb_info2[res_no2][2], Cb_info[res_no][2])

                                        #print(str(res_no) + "," + str(res_no2) + " phi = " + str(phi))
                                        omega = getDihedral(Ca_info[res_no][0], Cb_info[res_no][0], Cb_info2[res_no2][0], Ca_info2[res_no2][0],
                                                            Ca_info[res_no][1], Cb_info[res_no][1], Cb_info2[res_no2][1], Ca_info2[res_no2][1],
                                                            Ca_info[res_no][2], Cb_info[res_no][2], Cb_info2[res_no2][2], Ca_info2[res_no2][2]) 
                                        
                                        #print(str(res_no) + "," + str(res_no2) + " theta = " + str(theta))

                                        outFile.write(str(res_no) + ' ' + str(res_no2) + ' ' + str(phi) + ' ' + str(phi2) +  ' ' + str(theta) + ' ' + str(theta2) + ' ' + str(omega) + ' ' + chainIds[0] + ' ' + chainIds[1] + "\n")
                                except Exception as e:
                                        print('Error')
                outFile.close()       
                
def main():
    chainIds = [chain1, chain2]
    pdbList = os.listdir(decoyDir)
    for d in range(len(pdbList)):
        print(decoyDir + '/' + pdbList[d])
        pdbFile = decoyDir + '/' + pdbList[d]
        outFile = open(outPut + '/' + pdbList[d].split('.pdb')[0] + '.ori', 'w')
        cal_orientation(pdbFile, outFile, chainIds)
        
if __name__ == "__main__":
        main()


