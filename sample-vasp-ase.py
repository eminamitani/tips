from ase import Atoms
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.build import fcc111
from ase.calculators.vasp import Vasp

#read structure from CONTCAR
mol = read('CONTCAR', format='vasp')
#move to center
mol.translate([mol.cell[0][0]/2,mol.cell[1][1]/2,mol.cell[2][2]/2])
mol.wrap(pbc=[1,1,1])

#set supercell
#supercell size: 8x8 4-layer
x=7
y=7
z=4

vac=18.0
slab=fcc111('Au', a=4.9, size=(x,y,z), vacuum=vac)

adsorbates = Atoms(mol.symbols, mol.positions)
adsorbates.set_cell(slab.get_cell(), scale_atoms=False)
#set center to rotate the molecule
#for example Fe site
c=[]
for atom in adsorbates:
    if (atom.symbol == 'Fe'):
        c.append(atom.position)
#rotate 30 degree
angle=30.0
adsorbates.rotate(angle, 'z', center=c[0])

#find ontop-site
topLayers = []
top = max(slab.positions[:, 2])
bottom=min(slab.positions[:, 2])

for _ in slab.positions:
    if (abs(_[2] - top)) < 0.1:
        topLayers.append(_)

ontop = topLayers[0]
print("ontop:"+str(ontop))
#define adsorption sit
vector=ontop - c[0]
adsorbates.translate(vector)

#move upward 2A
dist=2.0
adsorbates.translate([0.0, 0.0, top+dist])

slab.extend(adsorbates)
slab.center(vacuum=vac, axis=2)

fix=FixAtoms(indices=[atom.index for atom in slab if abs(atom.position[2]-bottom) < 0.5])
slab.set_constraint(fix)

write("POSCAR.vasp", slab, sort=True, format='vasp', vasp5=True, direct=True)
calc = Vasp(xc='pbe',setups='recommended')
atoms=read('POSCAR.vasp',format='vasp')
calc.initialize(atoms)
calc.write_potcar()


