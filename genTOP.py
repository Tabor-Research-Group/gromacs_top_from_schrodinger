import argparse
import numpy as np
from ase.io import read

#test

def printGro(infile, outfile, atomSymbol):
    mol = read(infile)
    coord = mol.get_positions()

    gout = open(outfile, 'w')

    print('Gro File', file=gout)
    print(natoms, file=gout)
    for i in range(len(coord)):
        print("{0:>5d}{1:<5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}".format              (1, 'MOL', atomSymbol[i], i+1, coord[i][0]/10, coord[i][1]/10, coord[i][2]/10), file=gout)

    gout.close()



def printGro2(infile, outfile, atomSymbol, cell):
    mol = read(infile)
    coord = mol.get_positions()

    gout = open(outfile, 'w')

    print('Gro File', file=gout)
    print(natoms, file=gout)
    for i in range(len(coord)):
        print("{0:>5d}{1:<5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}".format              (1, 'MOL', atomSymbol[i], i+1, coord[i][0]/10, coord[i][1]/10, coord[i][2]/10), file=gout)
    print("{0:>10.5f}{1:>10.5f}{2:>10.5f}".format(cell[0][0]/10, cell[1][1]/10, cell[2][2]/10), file=gout)
    gout.close()



def grep_atom_types(filename, atomNum, masses):
    with open(filename, 'r') as f:
        lines = f.readlines()
        buffer = []
        log = False
        for line in lines:
            if line.startswith('OPLSAA'):
                buffer.append(line)
                log = True
            elif line.startswith(' Stretch'):
                buffer.append(line)
                log = False
            elif log:
                buffer.append(line)

        tp = buffer[4:-4]
        atomName = [i.split()[0] for i in tp]
        atomType = [i.split()[1] for i in tp]
        atomSymbol = [i.split()[3] for i in tp]
        charges = [i.split()[4] for i in tp]
        sigma = [i.split()[5] for i in tp]
        epsilon = [i.split()[6] for i in tp]
        index = np.arange(1, len(atomName)+1).tolist()

        atom_dict = dict(zip(atomName, index))
        nb_dict = {}
        for i in range(len(atomName)):
            d = {atomType[i]:[atomSymbol[i], atomNum[i], masses[i], sigma[i], epsilon[i]]}
            nb_dict.update(d)

    return atomName, atomType, atomSymbol, charges, sigma, epsilon, atom_dict, nb_dict



def parseFF(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        buffer1 = []
        buffer2 = []
        buffer3 = []
        buffer4 = []
        log = False
        for line in lines:
            if line.startswith(' Stretch'):
                buffer1.append(line)
                log = True
            elif line.startswith(' Bending'):
                buffer1.append(line)
                log = False
            elif log:
                buffer1.append(line)
        for line in lines:
            if line.startswith(' Bending'):
                buffer2.append(line)
                log = True
            elif line.startswith(' proper Torsion'):
                buffer2.append(line)
                log = False
            elif log:
                buffer2.append(line)
        for line in lines:
            if line.startswith(' proper Torsion'):
                buffer3.append(line)
                log = True
            elif line.startswith(' improper Torsion'):
                buffer3.append(line)
                log = False
            elif log:
                buffer3.append(line)
        for line in lines:
            if line.startswith(' improper Torsion'):
                buffer4.append(line)
                log = True
            elif log:
                buffer4.append(line)

        tp = buffer1[1:-2]
        bond = [i.split()[:4] for i in tp]
        tp = buffer2[1:-2]
        angle = [i.split()[:5] for i in tp]
        tp = buffer3[1:-2]
        torsion = [i.split()[:8] for i in tp]
        tp = buffer4[1:-2]
        im_torsion = [i.split()[:5] for i in tp]


        return bond, angle, torsion, im_torsion



def dihedral_conv(torsions):
    dihedrals = []
    for i in torsions:
        v1 = float(i[4])
        v2 = float(i[5])
        v3 = float(i[6])
        v4 = float(i[7])
        c0 = (v2 + 0.5*(v1+v3))*8.368
        c1 = 0.5*(-1*v1 + 3*v3)*8.368
        c2 = (-1*v2 + 4*v4)*8.368
        c3 = -2*v3*8.368
        c4 = -4*v4*8.368
        c5 = 0

        new_tor = i[:4] + [c0, c1, c2, c3, c4, c5]
        dihedrals.append(new_tor)

    return dihedrals


def printTOP(outfile, natoms, atomType, atomSymbol, charges, masses, bonds, angles, dihedrals, im_torsions, atom_dict):

    fout = open(outfile, 'w')

    # Begin Writing Topology File
    print('; Topology file for MOL', file=fout)
    print('', file=fout)
    print('[ moleculetype ]', file=fout)
    print('; name     nrexcl', file=fout)
    print('  ', 'MOL', '     3', file=fout)
    print('', file=fout)

    # Print Atoms
    print('[ atoms ]', file=fout)
    print(';   nr     type  resnr residue    atom     cgnr    charge     mass', file=fout)
    for i in range(natoms):
        print('{0:>6d}{1:>10s}{2:>6d}{3:>7s}{4:>10s}{5:>9d}{6:>10.3f}{7:>10.3f}'.format              (i+1, 'opls_'+atomType[i], 1, 'MOL', atomSymbol[i], i+1, float(charges[i]), masses[i]), file=fout)


    # Print Bonds
    print('', file=fout)
    print('[ bonds ]', file=fout)
    print(';   ai    aj  funct        c0        c1', file=fout)
    for i in range(len(bonds)):
        print('{0:>6d}{1:>6d}{2:>7d}{3:>10.5f}{4:>10.1f}'.format              (atom_dict[bonds[i][0]], atom_dict[bonds[i][1]], 1, float(bonds[i][3])/10, float(bonds[i][2])*4.184*200), file=fout)

    # Print Angles
    print('', file=fout)
    print('[ angles ]', file=fout)
    print(';   ai    aj    ak   funct    theta0       k0', file=fout)
    for i in range(len(angles)):
        print('{0:>6d}{1:>6d}{2:>6d}{3:>8d}{4:>10.3f}{5:>10.1f}'.format              (atom_dict[angles[i][0]], atom_dict[angles[i][1]], atom_dict[angles[i][2]], 1, float(angles[i][4])*8.368, float(angles[i][3])*4.184), file=fout)


    # Print Dihedrals
    print('', file=fout)
    print("[ dihedrals ]", file=fout)
    print(";   ai    aj    ak    al funct         c0           c1          c2          c3          c4          c5(kj/mol)", file=fout)
    for i in range(len(dihedrals)):
        print("{0:>6d}{1:>6d}{2:>6d}{3:>6d}{4:>6d}{5:>12.6f}{6:>12.6f}{7:>12.6f}{8:>12.6f}{9:>12.6f}{10:>12.6f}".format              (atom_dict[dihedrals[i][0]], atom_dict[dihedrals[i][1]], atom_dict[dihedrals[i][2]], atom_dict[dihedrals[i][3]], 3, dihedrals[i][4], dihedrals[i][5],               dihedrals[i][6], dihedrals[i][7], dihedrals[i][8], dihedrals[i][9]), file=fout)


    # Print Improper Dihedrals
    print('', file=fout)
    print("[ dihedrals ]", file=fout)
    print(";improper dihedrals", file=fout)
    print(";    ai   aj    ak    al funct    phi(deg)  k(kJ/mol)   multi", file=fout)
    for i in range(len(im_torsions)):
        print("{0:>6d}{1:>6d}{2:>6d}{3:>6d}{4:>6d}{5:>12.1f}{6:>11.5f}{7:>8d}".format              (atom_dict[im_torsions[i][0]], atom_dict[im_torsions[i][1]], atom_dict[im_torsions[i][2]], atom_dict[im_torsions[i][3]],               1, 180, float(im_torsions[i][4])*8.368, 2), file=fout)


    fout.close()



def printNB(nb_dict):
    fout = open('ffnonbonded.itp', 'w')
    print(' [ atomtypes ]', file=fout)
    print('; full atom descriptions are available in ffoplsaa.atp', file=fout)
    print('; name  bond_type       mass    charge ptype   sigma epsilon', file=fout)
    for i in nb_dict:
        print("{0:>8s}{1:>6s}{2:>6d}{3:>10.5f}{4:>10.3f}{5:>6s}{6:>8.4f}{7:>8.4f}".format              ('opls_'+i, nb_dict[i][0], nb_dict[i][1], nb_dict[i][2], 0.000, 'A', float(nb_dict[i][3])/10, float(nb_dict[i][-1])*4.184), file=fout)


    fout.close()



def printITP():
    fout = open('topol.top', 'w')
    print('; top', file=fout)
    print('', file=fout)
    print('[ defaults ]', file=fout)
    print('; nbfunc    comb-rule   gen-pairs   fudgeLJ   fudgeQQ', file=fout)
    print('       1            3         yes       0.5       0.5', file=fout)
    print('', file=fout)
    print('#include "ffnonbonded.itp"', file=fout)
    print('#include "mol.top"', file=fout)
    print('', file=fout)
    print('[system]', file=fout)
    print('System', file=fout)
    print('', file=fout)
    print('[molecules]', file=fout)
    print('', file=fout)
    print('MOL 1', file=fout)
    print('', file=fout)



if __name__ == '__main__':
    # Parse Command-line Input
    parser = argparse.ArgumentParser(description='Generate random crystals from input(.log, .gjf or .xyz) file')
    parser.add_argument('-f', nargs=1, help='Input geometry file (xyz, gro)', required=True)
    parser.add_argument('-l', nargs=1, help='Input log file', required=True)
    parser.add_argument('-o', nargs=1, help='Output gromacs gro file', required=True)
    parser.add_argument('-t', nargs=1, help='Topology file', required=True)
    args = parser.parse_args()

    if vars(args)['l'][0] != 'NULL' and vars(args)['f'][0] != 'NULL':
        geom_input = vars(args)['f'][0]
        log_input = vars(args)['l'][0]
        new_gro_output = vars(args)['o'][0]
        top_output = vars(args)['t'][0]

        bonds, angles, torsions, im_torsions = parseFF(log_input)

        mol = read(geom_input)
        masses = mol.get_masses()
        natoms = len(masses)
        atomNum = mol.get_atomic_numbers()
        atomName, atomType, atomSymbol, charges, sigma, epsilon, atom_dict, nb_dict = grep_atom_types(log_input, atomNum, masses)
        dihedrals = dihedral_conv(torsions)

        printNB(nb_dict)
        printITP()
        printTOP(top_output, natoms, atomType, atomSymbol, charges, masses, bonds, angles, dihedrals, im_torsions, atom_dict)

        if vars(args)['f'][0][-3:] == 'gro':
            cell = mol.get_cell()
            printGro2(geom_input, new_gro_output, atomSymbol, cell)
        else:
            printGro(geom_input, new_gro_output, atomSymbol)



