import numpy as np
import matplotlib.pyplot as plt

# Define beam parameters
L = 1.0         # Length of the beam (m)
w = 0.01        # Width of the beam (m)
h = 0.02        # Height of the beam (m)
P = 100.0       # Magnitude of the tip load (N)
E = 2.0e11      # Young's modulus of the beam material (Pa)
nu = 0.3        # Poisson's ratio of the beam material

# Define finite element parameters
n_elements = 10         # Number of finite elements
n_nodes = n_elements+1  # Number of nodes
dx = L/n_elements       # Length of each element

# Initialize arrays
K = np.zeros((n_nodes,n_nodes))  # Stiffness matrix
F = np.zeros(n_nodes)            # Load vector
u = np.zeros(n_nodes)            # Displacement vector

# Define element stiffness matrix
def kelement(E, nu, w, h, dx):
    # Calculate constants
    A = w*h
    I = (w*h**3)/12
    k1 = 12*E*I/(dx**3*(1-nu**2))
    k2 = 6*E*I/(dx**2*(1-nu**2))
    k3 = 4*E*I/(dx*(1-nu**2))
    k4 = 2*E*I*(1+nu)/(dx*(1-nu**2))
    k5 = A*E/(1-nu**2)
    
    # Assemble stiffness matrix
    K = np.array([[k1, k2, -k3, k1, k2, k3],
                  [k2, k1+2*k4, k5, k2, k1-k4, -k5],
                  [-k3, k5, 2*k1, k3, -k5, k1-k4],
                  [k1, k2, k3, k1, k2, -k3],
                  [k2, k1-k4, -k5, k2, k1+2*k4, k5],
                  [k3, -k5, k1-k4, -k3, k5, 2*k1]])
    return K

# Assemble global stiffness matrix
for i in range(n_elements):
    ke = kelement(E, nu, w, h, dx)
    K[i:i+6, i:i+6] += ke

# Apply tip load
F[-1] = P

# Solve system of equations
u = np.linalg.solve(K, F)

# Plot results
x = np.linspace(0, L, n_nodes)
plt.plot(x, u)
plt.xlabel('Position (m)')
plt.ylabel('Displacement (m)')
plt.title('Cantilever beam with tip load')
plt.show()
