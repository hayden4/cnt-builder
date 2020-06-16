import math
import sys

class CarbonNanotube:
    
    def __init__(self, s, n, L):
        self.s = s
        self.num_atoms = n+1
        self.L = L
        self.atoms = []
        self.bonds = []


    def buildCNT(self):
        dt = (2 * math.pi) / self.num_atoms
        radius = self.s * math.sqrt(3 / (2 * (1 - math.cos(dt))))
        for i in range(self.num_atoms):
            self.atoms.append([radius * math.cos(i * dt), radius * math.sin(i * dt), 0])


        z = 0
        while z < self.L:
            dt = (2 * math.pi) / self.num_atoms

            # sublayer 1
            z += self.s/2
            if z > self.L: break
            for i in range(self.num_atoms-1):
                self.atoms.append([radius * math.cos(i*dt + (dt/2)), radius * math.sin(i*dt + (dt/2)), z])
                self.bonds.append([len(self.atoms) - self.num_atoms    , len(self.atoms)])
                self.bonds.append([len(self.atoms) - self.num_atoms + 1, len(self.atoms)])
            self.atoms.append([radius * math.cos(-dt/2), radius * math.sin(-dt/2), z])
            self.bonds.append([len(self.atoms) - 2*self.num_atoms+1, len(self.atoms)])
            self.bonds.append([len(self.atoms) -   self.num_atoms  , len(self.atoms)])

            # sublayer 2
            z += self.s
            if z > self.L: break
            for i in range(self.num_atoms):
                self.atoms.append([radius * math.cos(i*dt + (dt/2)), radius * math.sin(i*dt + (dt/2)), z])
                self.bonds.append([len(self.atoms) - self.num_atoms, len(self.atoms)])

            # sublayer 3
            z += self.s/2
            if z > self.L: break
            self.atoms.append([radius, 0, z])
            self.bonds.append([len(self.atoms)-1             , len(self.atoms)])
            self.bonds.append([len(self.atoms)-self.num_atoms, len(self.atoms)])
            for i in range(1, self.num_atoms):
                self.atoms.append([radius * math.cos(i*dt), radius * math.sin(i*dt), z])
                self.bonds.append([len(self.atoms)-self.num_atoms-1, len(self.atoms)])
                self.bonds.append([len(self.atoms)-self.num_atoms  , len(self.atoms)])

            # sublayer 4
            z += self.s
            if z > self.L: break
            for i in range(self.num_atoms):
                self.atoms.append([radius * math.cos(i*dt), radius * math.sin(i*dt), z])
                self.bonds.append([len(self.atoms)-self.num_atoms, len(self.atoms)])

    def translate(self, r):
        for i in range(len(self.atoms)):
            self.atoms[i][0] += r


    def rotateCNT(self, angle):
        dy = self.L * math.sin(angle) / 2
        dz = self.L * (1 - math.cos(angle)) / 2
        for i in range(len(self.atoms)):
            pos = self.atoms[i]
            self.atoms[i] = [
                pos[0],
                pos[1] * math.cos(angle) - pos[2] * math.sin(angle) + dy,
                pos[1] * math.sin(angle) + pos[2] * math.cos(angle) + dz]




# builds a armchair (n,n) CNT
def main():
    separation = float(sys.argv[1])

    # CNT Configuration
    # s = C-C bond length
    # n = number of carbon atoms in ring -1
    # L = length of CNT

    radius = 1.0
    n = 10
    L = 25 * radius

    # Set bond length from parameters
    s = radius * math.sqrt(2*(1-math.cos(2*math.pi/(n+1)))/3)

    # Build CNT 1
    cnt1 = CarbonNanotube(s, n, L)
    cnt1.buildCNT()

    # Build CNT 2
    cnt2 = CarbonNanotube(s, n, L)
    cnt2.buildCNT()
    cnt2.translate(20)
    cnt2.rotateCNT(math.pi/2)

    # Combine tubes
    atoms = [a for a in cnt1.atoms]
    bonds = [b for b in cnt1.bonds]

    atomid_offset = len(atoms)
    for i in range(len(cnt2.atoms)):
        atoms.append(cnt2.atoms[i])

    for i in range(len(cnt2.bonds)):
        bonds.append([cnt2.bonds[i][0]+atomid_offset, cnt2.bonds[i][1]+atomid_offset])

    # Generate initial ring

    outfile = open('cnt.lammps', 'w')
    outfile.write("LAMMPS Description\n\n")
    outfile.write("\t" + str(len(atoms)) + " atoms\n")
    outfile.write("\t" + str(len(bonds)) + " bonds\n")
    outfile.write("\t1 atom types\n")
    outfile.write("\t1 bond types\n")
    outfile.write("\t{a} {b} xlo xhi\n".format(a=-1000, b=1000))
    outfile.write("\t{a} {b} ylo yhi\n".format(a=-1000, b=1000))
    outfile.write("\t{a} {b} zlo zhi\n".format(a=-1000, b=1000))
    outfile.write("\nMasses\n\n")
    outfile.write("1 1\n")
    outfile.write("\nAtoms\n\n")

    for i in range(len(atoms)):
        outfile.write("{i} 1 1 0 {x} {y} {z}\n".format(i=i+1, x=atoms[i][0], y=atoms[i][1], z=atoms[i][2]))

    outfile.write("\nBonds\n\n")

    for i in range(len(bonds)):
        outfile.write("{i} 1 {a} {b}\n".format(i=i+1, a=bonds[i][0], b=bonds[i][1]))

    outfile.close()
        



if __name__ == '__main__':
    main()
