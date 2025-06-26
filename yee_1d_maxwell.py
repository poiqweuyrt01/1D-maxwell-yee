# Yee 1D Maxwell's Equations Simulation



import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display, clear_output
import time

def yee_mesh(xl, xr, N): #creates the meshnodes for periodic boundary conditions
    #inputs:
    # xl: the left endpoint of the interval
    # xr: the right endpoint of the interval
    # N: the number of mesh nodes for the electric field grid
    #########################################################
    #outputs:
    # h: the spacing between meshnodes
    # x_E: the meshnodes for the electric field
    # x_H: the meshnodes for the magnetic field
    
    h = (xr-xl)/(N-1) 
    x_E = np.linspace(xl, xr, N)
    x_H = np.linspace(xl+h/2, xr+h/2, N)
    return h, x_E, x_H

def yee1D(x_E, x_H, h, E0, H0, t0, T, cfl): #calculate the electric and magnetic fields for a 
                                            #non-disapative, lossless, isotropic material
    #inputs:
    # x_E: the meshnodes for the electric field
    # x_H: the meshnodes for the magnetic field    
    # h: the spacing between meshnodes
    # E0: the initial condition for the electric field
    # H0: the initial condition for the magnetic field
    # t0: the initial time
    # T: the final time to compute the electric field
    # cfl: the mesh ratio
    
    E = np.copy(E0)
    H = np.copy(H0)
    t = t0
    while t < T:
        H_bdy = H[0] - cfl*(E[1] - E[0])
        H[1:-1] = H[1:-1] - cfl*(E[2:] - E[1:-1])
        H[0] = H_bdy
        H[-1] = H_bdy
        E_bdy = E[0] - cfl*(H[0] - H[-2])
        E[1:-1] = E[1:-1] - cfl*(H[1:-1] - H[:-2])
        E[0] = E_bdy
        E[-1] = E_bdy
        t += cfl*h
    return H, E, t

cfl = 0.9
E0 = lambda x: np.cos(x)  ##it is the same thing as def E0(x): return np.cos(x)
H0 = lambda x, cfl, h: np.sin(x)*np.sin(-cfl*h/2)
E_ex = lambda x, t: np.cos(x)*np.cos(t)
H_ex = lambda x, t: np.sin(x)*np.sin(t)
h, x_E, x_H = yee_mesh(0, 2*np.pi, 100)
H, E, t = yee1D(x_E, x_H, h, E0(x_E), H0(x_H,cfl,h), 0, 10, cfl)

plt.plot(x_E, E, color='red', marker='o', label='Computed Electric Field')
plt.plot(x_E, E_ex(x_E, t), color='blue', label='Exact Electric Field')
plt.plot(x_H, H, color='orange', marker='o', label='Computed Magnetic Field')
plt.plot(x_H, H_ex(x_H, t-cfl*h/2), color='green', label='Exact Magnetic Field')
plt.xlabel('$x$')
plt.ylabel('Field Strength')
plt.legend()
plt.title('Electric and Magnetic Field at $t=10$')

N = [2**k+1 for k in range(5, 15)]
N = np.array(N)
h_vals = 2*np.pi/(N-1)
cfl = 0.9
e = np.zeros((len(N),2))
E_ex = lambda x, t: np.cos(x)*np.cos(t)
H_ex = lambda x, t: np.sin(x)*np.sin(t)
for i in range(10):
    h, x_E, x_H = yee_mesh(0, 2*np.pi, N[i])
    E0 = lambda x: np.cos(x)
    H0 = lambda x, cfl, h: np.sin(x)*np.sin(-cfl*h/2)
    H, E, t = yee1D(x_E, x_H, h, E0(x_E), H0(x_H,cfl,h), 0, 10, cfl)
    E_exact = E_ex(x_E, t)
    H_exact = H_ex(x_H, t-cfl*h/2)
    e[i,0] = np.linalg.norm(E-E_exact, np.inf)
    e[i,1] = np.linalg.norm(H-H_exact, np.inf)
plt.loglog(h_vals, e[:,0], color='red', marker='o', label=r'$||E-E^*||_{\infty}$')
plt.loglog(h_vals, e[:,1], color='orange', marker='o', label=r'$||H-H^*||_{\infty}$')
plt.loglog(h_vals, h_vals**2, color='blue', label='$h^2$')
plt.legend()
plt.title("Error versus $h$")
plt.xlabel('$h$')
plt.ylabel('Error')

def yee_mesh_np(xl, xr, N): #creates the yee mesh for non-periodic boundary conditions
    #inputs:
    # xl: the left endpoint of the interval
    # xr: the right endpoint of the interval
    # N: the number of mesh nodes for the electric field grid
    #########################################################
    #outputs:
    # h: the spacing between meshnodes
    # x_E: the meshnodes for the electric field
    # x_H: the meshnodes for the magnetic field
    
    h = (xr-xl)/ (N-1)
    x_E = np.linspace(xl, xr, N)
    x_H = np.linspace(xl+h/2, xr-h/2, N-1)
    return h, x_E, x_H

L = 200
N_steps = 800
h, x_E, x_H = yee_mesh_np(0, L, 200)
ey = np.zeros((len(x_E), N_steps + 1), dtype=np.float64)
hz = np.zeros((len(x_H), N_steps + 1), dtype=np.float64)
d = L / 2
tau = 12
cfl = 0.9
a = 1

fig = plt.figure(figsize=(10, 6))

ey[:, 0] = a * np.exp(-0.5 * ((x_E - d) / tau) ** 2) 
# a is the height 
# d is the position of the center of the peak (the mean of the Gaussian distribution)
#tau is the standard deviation of the Gaussian distribution, controlling the width of the curve.

for n in range(N_steps):
    hz[:, n] = hz[:, n - 1]  - cfl * (ey[1:, n] - ey[:-1, n])
    ey[1:-1, n + 1] = ey[1:-1, n] - cfl * (hz[1:, n] - hz[:-1, n])
    plt.clf()
    plt.plot(x_E, ey[:, n], label='Electric Field')
    plt.plot(x_H, hz[:, n], label='Magnetic Field')
    plt.ylim([-1, 1])
    plt.title('Electric and Magnetic Fields with R')
    plt.legend()
    display(fig)
    clear_output(wait=True)
    time.sleep(0.01)
plt.close()

L = 200
N_steps = 800
h, x_E, x_H = yee_mesh_np(0, L, 201)
ey = np.zeros((len(x_E), N_steps + 1), dtype=np.float64)
hz = np.zeros((len(x_H), N_steps + 1), dtype=np.float64)
d = L / 2
tau = 12
cfl = 0.9

fig = plt.figure(figsize=(10, 6))

ey[:, 0] = np.exp(-0.5 * ((x_E - d) / tau) ** 2)

for n in range(N_steps):
    hz[:, n] = hz[:, n - 1]  - cfl * (ey[1:, n] - ey[:-1, n])
    ey[1:-1, n + 1] = ey[1:-1, n] - cfl * (hz[1:, n] - hz[:-1, n])
    ey[-1, n + 1] = ey[-2, n] + ((cfl - 1)/(cfl + 1)) * (ey[-2, n + 1] - ey[-1, n])
    plt.clf()
    plt.plot(x_E, ey[:, n], label='Electric Field')
    plt.plot(x_H, hz[:, n], label='Magnetic Field')
    plt.ylim([-1, 1])
    plt.title('Electric and Magnetic Fields with ABC')
    plt.legend()
    display(fig)
    clear_output(wait=True)
    time.sleep(0.01)
plt.close()

L = 200
N_steps = 800
h, x_E, x_H = yee_mesh_np(0, L, 201)
ey = np.zeros((len(x_E), N_steps + 1), dtype=np.float64)
hz = np.zeros((len(x_H), N_steps + 1), dtype=np.float64)
d = L / 2
tau = 12
c = 1

varepsilon = 1
mu = 1 / (varepsilon * c ** 2)
sigma = 0.005
dt = 1
dx = 1

fig = plt.figure(figsize=(10, 6))

ey[:, 0] = np.exp(-0.5 * ((x_E - d) / tau) ** 2)

for n in range(N_steps):
    hz[:, n] = hz[:, n - 1] -  dt / (mu * dx) * (ey[1:, n] - ey[:-1, n])
    ey[1:-1, n + 1] = (2 * varepsilon - sigma * dt) / (2 * varepsilon + sigma * dt) * ey[1:-1, n] \
                    - (2 * dt / ((2 * mu + sigma * dt) * dx))  * (hz[1:, n] - hz[:-1, n])
    plt.clf()
    plt.plot(x_E, ey[:, n], label='Electric Field')
    plt.plot(x_H, hz[:, n], label='Magnetic Field')
    plt.ylim([-1, 1])
    plt.title('Electric and Magnetic Fields with lossy media')
    plt.legend()
    display(fig)
    clear_output(wait=True)
    time.sleep(0.01)
plt.close()

#use Newton's method to find the solution to -1.5tan(w) = tan(1.5w) between 3pi/2 and 5pi/3
x = 0
x1 = 19*np.pi/12
while abs(x1-x)>1e-12:
    x = x1
    x1 = x - (np.tan(1.5*x)+1.5*np.tan(x))/(1.5/(np.cos(1.5*x)**2)+1.5/(np.cos(x)**2))
x1

N_steps = 800
h, x_E, x_H = yee_mesh_np(-1, 1, 201)
ey = np.zeros((len(x_E), N_steps + 1), dtype=np.float64)
hz = np.zeros((len(x_H), N_steps + 1), dtype=np.float64)

def eps(x, eps1, eps2):
    #inputs:
    # x: a vector containing meshnodes
    # eps1: the permittivity for x<=0
    # eps2: the permittivity for x>0
    ##################################
    #outputs:
    # eps: a vector containing the value of permittivity at each node
    eps = np.where(x<=0, eps1, eps2)
    return eps

ep = eps(x_E, 1, 2.25)
cfl = 0.9
mu = 1
n1 = 1
n2 = 1.5
w = 5.0721811618
A1 = (n2*np.cos(n2*w))/(n1*np.cos(n1*w))
B1 = A1*np.exp(-2*n1*w*1j)
A2 = np.exp(-w*(n1+n2)*1j)
B2 = A2*np.exp(2*n2*w*1j)

def e_exact(x, t, n1, n2, A1, B1, A2, B2):
    e1 = -(A1*np.exp(n1*w*x*1j)-B1*np.exp(-n1*w*x*1j))*np.exp(w*t*1j)
    e2 = -(A2*np.exp(n2*w*x*1j)-B2*np.exp(-n2*w*x*1j))*np.exp(w*t*1j)
    e_ex = np.where(x<=0, e1, e2)
    return e_ex

def h_exact(x, t, n1, n2, A1, B1, A2, B2):
    h1 = n1*(A1*np.exp(n1*w*x*1j)+B1*np.exp(-n1*w*x*1j))*np.exp(w*t*1j)
    h2 = n2*(A2*np.exp(n2*w*x*1j)+B2*np.exp(-n2*w*x*1j))*np.exp(w*t*1j)
    h_ex = np.where(x<=0, h1, h2)
    return h_ex

ey[:,0] = e_exact(x_E, 0, n1, n2, A1, B1, A2, B2)
hz[:,0] = h_exact(x_H, -cfl*h/2, n1, n2, A1, B1, A2, B2)

fig = plt.figure(figsize=(10, 6))

for n in range(1,N_steps):
    hz[:,n] = hz[:,n-1] - 1/mu * cfl * (ey[1:,n-1]-ey[:-1,n-1])
    ey[1:-1,n] = ey[1:-1,n-1] - 1/ep[1:-1] * cfl * (hz[1:,n]-hz[:-1,n])
    plt.clf()
    plt.plot(x_E, ey[:, n], label='Electric Field')
    plt.plot(x_H, hz[:, n], label='Magnetic Field')
    plt.ylim([-3, 3])
    plt.title('Electric and Magnetic Fields')
    plt.legend()
    display(fig)
    clear_output(wait=True)
    time.sleep(0.01)
plt.close()

def yee1D_varicoeff(x_E, x_H, h, eps, mu, E0, H0, t0, T, cfl):
    #inputs:
    # x_E: the meshnodes for the electric field
    # x_H: the meshnodes for the magnetic field
    # h: the spacing between meshnodes
    # eps: a vector containing the value of permittivity at each electric field node
    # mu: a scalar, the value of the magnetic permiability
    # E0: the initial condition for the electric field
    # H0: the initial condition for the magnetic field
    # t0: the initial time
    # T: the final time to compute the electric field
    # cfl: the mesh ratio
    ################################################################################
    #outputs:
    # Hz: the computed magnetic field at time t - 1/2 dt
    # Ey: the computed electric field at time t
    # t: the computed final time
    
    Ey = np.copy(E0)
    Hz = np.copy(H0)
    t = t0
    while t < T:
        Hz = Hz - 1/mu*cfl*(Ey[1:] - Ey[:-1])
        Ey[1:-1] = Ey[1:-1] - 1/eps[1:-1]*cfl*(Hz[1:] - Hz[:-1])
        t += cfl*h
    return Hz, Ey, t

N = np.array([2**k+1 for k in range(5, 15)])
h_vals = 2*np.pi/(N-1)
cfl = 10/13
mu = 1
n1 = 1
n2 = 1.5
w = 5.0721811618
A1 = (n2*np.cos(n2*w))/(n1*np.cos(n1*w))
B1 = A1*np.exp(-2*n1*w*1j)
A2 = np.exp(-w*(n1+n2)*1j)
B2 = A2*np.exp(2*n2*w*1j)
e = np.zeros((len(N),2))
for i in range(10):
    h, x_E, x_H = yee_mesh_np(-1, 1, N[i])
    ep = eps(x_E, 1, 2.25)
    Ey0 = e_exact(x_E, 0, n1, n2, A1, B1, A2, B2)
    Hz0 = h_exact(x_H, -cfl*h/2, n1, n2, A1, B1, A2, B2)
    H, E, t = yee1D_varicoeff(x_E, x_H, h, ep, 1, Ey0, Hz0, 0, 2*np.pi, cfl)
    E_exact = e_exact(x_E, t, n1, n2, A1, B1, A2, B2)
    H_exact = h_exact(x_H, t-cfl*h/2, n1, n2, A1, B1, A2, B2)
    e[i,0] = np.linalg.norm(E-E_exact, 2)*np.sqrt(h_vals[i])
    e[i,1] = np.linalg.norm(H-H_exact, 2)*np.sqrt(h_vals[i])
plt.loglog(h_vals, e[:,0], color='red', marker='o', label=r'$||E-E^*||_{l^2}$')
plt.loglog(h_vals, e[:,1], color='orange', marker='o', label=r'$||H-H^*||_{l^2}$')
plt.loglog(h_vals, h_vals, color='blue', label='$h$')
plt.legend()
plt.title("Error versus $h$")
plt.xlabel('$h$')
plt.ylabel('Error')

