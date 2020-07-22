import torch
import numpy as np

#cutoff function, return scalar
def fc(R, Rc):
    return torch.where(R<Rc, torch.cos(np.pi*R/Rc)*0.5+1, torch.tensor(0.0))

def fctahn(R, Rc):
    return torch.tanh(1-R/Rc)**3

def g5(xi,yi,zi,xj,yj,zj, xk,yk,zk,Rc, lamb,eta, zeta):
    Rij=torch.sqrt(((xi-xj)**2+(yi-yj)**2+(zi-zj)**2))
    Rik = torch.sqrt(((xi - xk) ** 2 + (yi - yk) ** 2 + (zi - zk) ** 2))
    cosijk = ((xi - xj) * (xi - xk) + (yi - yj) * (yi - yk) + (zi - zj) * (zi - zk)) / Rij / Rik
    factor=torch.pow(2.0,1-zeta)*torch.exp(-eta*(Rij*Rij+Rik*Rik))*torch.pow(lamb*cosijk+1.0, zeta)
    return factor*fc(Rij,Rc)*fc(Rik,Rc)

def g5_deriv_xi(xi,yi,zi,xj,yj,zj, xk,yk,zk,Rc, lamb,eta, zeta):
    deriv=0.0
    Rij=torch.sqrt(((xi-xj)**2+(yi-yj)**2+(zi-zj)**2))
    Rik = torch.sqrt(((xi - xk) ** 2 + (yi - yk) ** 2 + (zi - zk) ** 2))
    cosijk = ((xi - xj) * (xi - xk) + (yi - yj) * (yi - yk) + (zi - zj) * (zi - zk)) / Rij / Rik
    exp=torch.exp(-eta*(Rij*Rij+Rik*Rik))
    pow=torch.pow(1.0+lamb*cosijk,zeta)
    inner=(xi - xj) * (xi - xk) + (yi - yj) * (yi - yk) + (zi - zj) * (zi - zk)
    dcosijk=1.0/Rij/Rik*((xi-xj)*(1.0-inner/Rij/Rij)+(xi-xk)*(1.0-inner/Rik/Rik))
    scale=torch.pow(2.0,1.0-zeta)
    deriv+=lamb*zeta*torch.pow(1.0+lamb*cosijk,zeta-1.0)*dcosijk*exp*fc(Rij,Rc)*fc(Rik,Rc)*scale

    dexp=-2.0*eta*((xi-xj)+(xi-xk))*torch.exp(-eta*(Rij*Rij+Rik*Rik))
    deriv+=pow*dexp*fc(Rij,Rc)*fc(Rik,Rc)*scale

    dfcRij=-0.5*np.pi/Rc*(xi-xj)/Rij*torch.sin(np.pi*Rij/Rc)
    deriv+=pow*exp*dfcRij*fc(Rik,Rc)*scale

    dfcRik = -0.5 * np.pi / Rc * (xi - xk) / Rik * torch.sin(np.pi * Rik / Rc)
    deriv+=pow*exp*fc(Rij,Rc)*dfcRik*scale

    return deriv

#just test
xi=torch.tensor([0.0,1.0],requires_grad=True)
yi=torch.tensor([0.0,1.5],requires_grad=True)
zi=torch.tensor([0.0,1.2],requires_grad=True)
xj=torch.tensor([1.0,2.1],requires_grad=True)
yj=torch.tensor([1.0,2.2],requires_grad=True)
zj=torch.tensor([1.0,2.5],requires_grad=True)
xk=torch.tensor([1.2,2.2],requires_grad=True)
yk=torch.tensor([1.4,2.6],requires_grad=True)
zk=torch.tensor([1.5,2.5],requires_grad=True)
lamb=torch.tensor(1.0)
zeta=torch.tensor(1.0)
Rc=torch.tensor(6.0)
Rs=torch.tensor(0.0)
eta=torch.tensor(0.0)
g5=torch.sum(g5(xi,yi,zi,xj,yj,zj, xk,yk,zk,Rc, lamb, eta,zeta))
g5.backward()
print("derivative of g5")
print("xi grad")
print(xi.grad)
print("direct differentiation")
print(g5_deriv_xi(xi,yi,zi,xj,yj,zj, xk,yk,zk,Rc, lamb,eta, zeta))

print(yi.grad)
print(zi.grad)
print(xj.grad)
print(yj.grad)
print(zj.grad)
print(xk.grad)
print(yk.grad)
print(zk.grad)