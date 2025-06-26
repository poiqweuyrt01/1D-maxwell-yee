# Yee Algorithm for 1D Maxwell’s Equations

## 📚 Background and Motivation

This project explores the one-dimensional Maxwell’s equations for electromagnetic wave propagation in media that vary only in the \( x \)-direction. The governing equations under a linear medium and transverse electric (TE) mode are:

$$
\mu \frac{\partial H_z}{\partial t} + \frac{\partial E_y}{\partial x} = 0, \quad
\frac{\partial D}{\partial t} + \frac{\partial H_z}{\partial x} = 0
$$

Where:

- \( H_z \): z-component of the magnetic field  
- \( E_y \): y-component of the electric field  
- \( D = \varepsilon E_y + P \): Electric displacement field with polarization \( P \)

We consider the following boundary conditions:

- **Periodic**: \( x \in [0, 2\pi] \)
- **Neumann** and **Absorbing (ABC)**: \( x \in [0, 200] \)

---

## ⚙️ Numerical Method: Yee Algorithm (FDTD)

The Finite-Difference Time-Domain (FDTD) method is implemented via the Yee scheme to solve the 1D Maxwell equations.

### Reduced 1D Maxwell’s Equations (TE Mode):

$$
\frac{\partial E_y}{\partial t} = -\frac{1}{\varepsilon} \frac{\partial H_z}{\partial x}, \quad
\frac{\partial H_z}{\partial t} = -\frac{1}{\mu} \frac{\partial E_y}{\partial x}
$$

### Discretized Form (Yee Scheme):

$$
E_i^{k+1} = E_i^k - \frac{c \Delta t}{\Delta x} (H_{i+1/2}^{k+1/2} - H_{i-1/2}^{k+1/2})
$$

$$
H_{i+1/2}^{k+1/2} = H_{i+1/2}^{k-1/2} - \frac{c \Delta t}{\Delta x} (E_{i+1}^k - E_i^k)
$$

Where \( c = \frac{1}{\sqrt{\mu \varepsilon}} \) is the wave speed.

---

## 📐 Consistency & Stability

### 🔎 Consistency

- Local truncation error is \( \mathcal{O}(\Delta t^2) \)
- Method shows second-order convergence

### 🧮 Stability

- Stable under CFL condition: \( \Delta t < \Delta x \)
- Unstable for "magic timestep" \( \Delta t = \Delta x \)

---

## 🧪 Simulation Results

- **Standing Mode**: Numerical solutions converge to analytical standing wave solutions
- **Dirichlet BCs**: Cause wave reflection at boundaries
- **ABC**: Absorbing boundary successfully minimizes reflections

---

## ✅ Advantages

- Second-order accuracy
- Mimics physical wave propagation
- Simple, explicit time-marching scheme

## ❌ Limitations

- Cannot impose Dirichlet boundary conditions simultaneously on both \( E \) and \( H \)
- Stability requires small time steps (CFL-limited)

---

## 📎 References

1. Schneider, J.B. *Understanding the Finite-Difference Time-Domain Method*, 2010  
   [http://www.eecs.wsu.edu/~schneidj/ufdtd](http://www.eecs.wsu.edu/~schneidj/ufdtd)

2. Taflove, A., Hagness, S.C., & Piket-May, M. *Computational Electromagnetics: The Finite-Difference Time-Domain Method*, Elsevier, 2005