import torch
import numpy as np
import ase.io
import ase.neighborlist

#cutoff function, using tanh^3 type
def fctahn(R, Rc):
    return torch.tanh(1-R/Rc)**3

def fc(R, Rc):
    return torch.where(R<Rc, torch.cos(np.pi*R/Rc)*0.5+1, torch.tensor(0.0,dtype=torch.double))

#type-2 symmetry function, return tensor when input is tensor
def g2(xi,yi,zi, xj,yj,zj, Rc, Rs, eta):
    R=torch.sqrt(((xi-xj)**2+(yi-yj)**2+(zi-zj)**2))
    return torch.exp(-eta*(R-Rs)**2) * fc(R,Rc)


Rc=torch.tensor(6.0)
Rs=torch.tensor(0.0)
eta=torch.tensor(0.0)

structure=ase.io.read('POSCAR',format='vasp')
#cutoff must be scalar value
#using Rc.item for this purpose
#ASE neighbor list does not include self,
#but count the atoms in image cell within cutoff
i_list, j_list, D_list = ase.neighborlist.neighbor_list(
    'ijD', structure, Rc.item())
print(i_list)
print(j_list)
print(D_list)

pos=structure.get_positions()
natom=len(pos)
#center
xi=torch.tensor([pos[i][0] for i in i_list],requires_grad=True,dtype=torch.double)
yi=torch.tensor([pos[i][1] for i in i_list],requires_grad=True,dtype=torch.double)
zi=torch.tensor([pos[i][1] for i in i_list],requires_grad=True,dtype=torch.double)
#neighbour
xj=torch.tensor([pos[j][0] for j in j_list],requires_grad=True,dtype=torch.double)
yj=torch.tensor([pos[j][1] for j in j_list],requires_grad=True,dtype=torch.double)
zj=torch.tensor([pos[j][1] for j in j_list],requires_grad=True,dtype=torch.double)

print(xi.shape)
g2array=torch.zeros([natom,1])
g2array=g2(xi,yi,zi, xj,yj,zj, Rc, Rs, eta)
print(g2array.shape)
print(np.count_nonzero(i_list ==0))

#make some scalar to obtain grad
g2sum=torch.sum(g2array)
print(g2sum)
g2sum.backward()
print(xi.grad)
print(xj.grad)


