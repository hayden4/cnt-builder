import math
import sys

class FCC_Tube:
    
    def __init__(self, radius, length, a):
        self.radius = radius
        self.length = length
        self.a = a
        self.l = 4*(int(radius/a)+1)
        self.m = 4*(int(radius/a)+1)
        self.n = 4*(int(length))
        self.atoms = []

    def buildFCC(self):
        for k in range(self.n):
            self.buildLayer(k)
        self.translate([-self.a*(self.l-1)/2, -self.a*(self.m-1)/2, 0])
        self.cutTube()


    def buildLayer(self, k):
        # Construct layer 1
        for j in range(self.m):
            # Primitive cubic row
            for i in range(self.l):
                self.atoms.append([self.a*i, self.a*j, self.a*k])
                if i < self.l-1 and j < self.m-1:
                    self.atoms.append([self.a*(i+1/2), self.a*(j+1/2), self.a*k])

        # Construct layer 2
        if k < self.n-1:
            for j in range(self.m):
                for i in range(self.l):
                    if i < self.l-1:
                        self.atoms.append([self.a*(i+1/2), self.a*j, self.a*(k+1/2)])
                    if j < self.m-1:
                        self.atoms.append([self.a*i, self.a*(j+1/2), self.a*(k+1/2)])


    def translate(self, vec):
        for i in range(len(self.atoms)):
            self.atoms[i] = [self.atoms[i][j] + vec[j] for j in range(3)]

    def rotate(self, angle):
        dy = self.length * math.sin(angle)/2
        dz = self.length * (1 - math.cos(angle))/2
        for i in range(len(self.atoms)):
            pos = self.atoms[i]
            self.atoms[i] = [
                pos[0],
                pos[1] * math.cos(angle) - pos[2] * math.sin(angle) + dy,
                pos[1] * math.sin(angle) + pos[2] * math.cos(angle) + dz]

    def cutTube(self):
        keep_atoms = []
        rsqd = pow(self.radius,2)
        for atom in self.atoms:
            if pow(atom[0],2) + pow(atom[1],2) < rsqd and atom[2] < self.length:
                keep_atoms.append(atom)
        self.atoms = keep_atoms
            




# Builds a FCC Crystal
def fcc_system(filename, separation, angle):
    radius = 1.0
    length = 25 * radius
    bond_length = 0.3

    fcc1 = FCC_Tube(radius, length, bond_length)
    fcc2 = FCC_Tube(radius, length, bond_length)

    fcc1.buildFCC()
    fcc2.buildFCC()

    fcc2.translate([separation, 0, 0])
    fcc2.rotate(angle)


    atoms = [a for a in fcc1.atoms]
    atomid_offset = len(atoms)
    for atom in fcc2.atoms:
        atoms.append(atom)

    # Write to file
    outfile = open(filename, 'w')
    outfile.write("LAMMPS Description\n\n")
    outfile.write("\t" + str(len(atoms)) + " atoms\n")
    outfile.write("\t1 atom types\n")
    outfile.write("\t{a} {b} xlo xhi\n".format(a=-1000, b=1000))
    outfile.write("\t{a} {b} ylo yhi\n".format(a=-1000, b=1000))
    outfile.write("\t{a} {b} zlo zhi\n".format(a=-1000, b=1000))
    outfile.write("\nMasses\n\n")
    outfile.write("1 1\n")
    outfile.write("\nAtoms\n\n")

    for i in range(atomid_offset):
        outfile.write("{i} 1 1 0 {x} {y} {z}\n".format(i=i+1, x=atoms[i][0], y=atoms[i][1], z=atoms[i][2]))
    
    for i in range(atomid_offset, len(atoms)):
        outfile.write("{i} 2 1 0 {x} {y} {z}\n".format(i=i+1, x=atoms[i][0], y=atoms[i][1], z=atoms[i][2]))

    outfile.close()
