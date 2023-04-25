'''
Finite Method on a Cantilever: ChatGPT edition
Author: Anupsa Swain, using ChatGPT
'''

# Prompt #1: write a python code for constructing a cantilever bar under weight at the end
# ChatGPT generated the following code:

import matplotlib.pyplot as plt
import numpy as np

# Define bar properties
L = 10.0  # Length of bar in meters
A = 0.1  # Cross-sectional area in square meters
E = 2e11  # Young's modulus in Pascals

# Define weight and position
W = 1.0e20  # Weight in Newtons
x_w = L  # Position of weight in meters

# Define number of elements
n = 100

# Define element length
dx = L / n

# Define element stiffness
k = E * A / dx

# Initialize arrays
u = np.zeros(n+1)
f = np.zeros(n+1)

# Set boundary condition
u[0] = 0.0

# Apply load at the end of the bar
f[-1] = W

# Construct stiffness matrix and load vector
K = np.zeros((n+1, n+1))
F = np.zeros(n+1)

for i in range(n):
    K[i,i] += k
    K[i+1,i+1] += k
    K[i,i+1] += -k
    K[i+1,i] += -k
    F[i] += f[i]
    F[i+1] += f[i+1]

print(K)
print(F)

# Solve system of equations
u = 10000*np.linalg.solve(K, F)

print(u)

# Plot deflection
x = np.linspace(0, L, n+1)
plt.plot(x, u)
plt.xlabel('Position along bar (m)')
plt.ylabel('Deflection (m)')
plt.title('Deflection of cantilever bar with weight at the end')
plt.show()
