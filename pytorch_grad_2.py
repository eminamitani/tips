import torch
import numpy as np

#just test
xi=torch.tensor([0.0,1.0],requires_grad=True)
yi=torch.tensor([0.0,1.5],requires_grad=True)
zi=torch.tensor([0.0,1.2],requires_grad=True)
xj=torch.tensor([1.0,2.1],requires_grad=True)
yj=torch.tensor([1.0,2.2],requires_grad=True)
zj=torch.tensor([1.0,2.5],requires_grad=True)
Rc=torch.tensor(6.0)
Rs=torch.tensor(0.0)
eta=torch.tensor(0.0)
print(xi)
#cutoff function, return scalar
def fc(R, Rc):
    return torch.where(R<Rc, torch.cos(np.pi*R/Rc)*0.5+1, torch.tensor(0.0))

def fctahn(R, Rc):
    return torch.tanh(1-R/Rc)**3


#type-2 symmetry function, return tensor when input is tensor
def g2(xi,yi,zi, xj,yj,zj, Rc, Rs, eta):
    R=torch.sqrt(((xi-xj)**2+(yi-yj)**2+(zi-zj)**2))
    return torch.exp(-eta*(R-Rs)**2) * fc(R,Rc)

#in order to obtain derivative for vector input (various atomic position etc)
#need to obtain scolar output to use autodifferentiation with grad
g2=torch.sum(g2(xi,yi,zi, xj,yj,zj, Rc, Rs, eta))
g2.backward()
print(xi.grad)
print(yi.grad)
print(zi.grad)
print(xj.grad)
print(yj.grad)
print(zj.grad)
