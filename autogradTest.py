#

import autograd.numpy as np
from autograd import grad

def g2(ri, rj, Rc, Rs, eta):
    R=np.sqrt(((ri[0]-rj[0])**2+(ri[1]-rj[1])**2+(ri[2]-rj[2])**2))
    fc = np.tanh(1.0 - R / Rc) ** 3

    return np.exp(-eta*(R-Rs)**2) * fc

grad(g2)
