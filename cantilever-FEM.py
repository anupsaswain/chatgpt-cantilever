import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

# Define the properties of the cantilever
L = 1.0     # Length of cantilever
E = 1.0e7   # Young's modulus
I = 1.0e-4  # Moment of inertia
w = 100.0   # Tip load

# A's added definition of analytical solution for comparison
def cantilytical(x):
    return w*(x**3)/(6*E*I) - w*L*x*x/(2*E*I)

# Define the boundary conditions
x0 = 0.0
y0 = 0.0
dydx0 = 0.0

# Define the domain of the solution
a = 0.0
b = L
n = 101   # Number of grid points
h = (b - a) / (n - 1)   # Grid spacing

# Define the element properties
ne = n - 1   # Number of elements
Le = L / ne  # Length of each element
ke = E * I / Le**3 * np.array([[12, 6*Le, -12, 6*Le],
                               [6*Le, 4*Le**2, -6*Le, 2*Le**2],
                               [-12, -6*Le, 12, -6*Le],
                               [6*Le, 2*Le**2, -6*Le, 4*Le**2]])

# Initialize the global stiffness matrix and load vector
K = lil_matrix((n, n))
F = np.zeros(n)

# Assemble the global stiffness matrix and load vector
for i in range(ne):
    K[i:i+4, i:i+4] += ke
    F[i:i+2] += w/2 * np.array([Le/2, Le/2])

# Apply boundary conditions
K = K[1:, 1:]
F = F[1:]
F[0] -= y0 * K[0, 0]
K = K[1:, 1:]
F = F - K[:, 0] * dydx0

# Solve for the deflection
y = np.zeros(n-1)
y[1:] = spsolve(K, F)

# Plot the solution
x = np.linspace(a, b, n)
plt.plot(x, np.concatenate(([y0], y))*1000, label='Deflection')
plt.xlabel('x (m)')
plt.ylabel('y (mm)')
plt.legend()
plt.show()
