import torch
import numpy as np

#just test
xi=torch.tensor(0.0,requires_grad=True)
yi=torch.tensor(0.0,requires_grad=True)
zi=torch.tensor(0.0,requires_grad=True)
xj=torch.tensor(1.0,requires_grad=True)
yj=torch.tensor(1.0,requires_grad=True)
zj=torch.tensor(1.0,requires_grad=True)
Rc=torch.tensor(6.0)
Rs=torch.tensor(0.0)
eta=torch.tensor(0.0)

#cutoff function, return scalar
def fc(R, Rc):
    return torch.where(R<Rc, torch.cos(np.pi*R/Rc)*0.5+1, torch.tensor(0.0))

def g2(xi,yi,zi, xj,yj,zj, Rc, Rs, eta):
    R=torch.sqrt(((xi-xj)**2+(yi-yj)**2+(zi-zj)**2))
    return torch.exp(-eta*(R-Rs)**2) * fc(R,Rc)

g2=g2(xi,yi,zi, xj,yj,zj, Rc, Rs, eta)
g2.backward()
print(xi.grad)
print(yi.grad)
print(zi.grad)
print(xj.grad)
print(yj.grad)
print(zj.grad)

'''
direct definition test

R=torch.sqrt(((xi-xj)**2+(yi-yj)**2+(zi-zj)**2))
fc=0.0
if R<Rc:
    fc=torch.cos(np.pi*R/Rc)*0.5+1
g2=torch.exp(-eta*(R-Rs)**2) * fc
g2.backward()
'''









