import numpy as np
import matplotlib.pyplot as plt

# Define the properties of the cantilever
L = 1.0     # Length of cantilever
E = 1.0e7   # Young's modulus
I = 1.0e-4  # Moment of inertia
w = 100.0   # Tip load

# Define the boundary conditions
x0 = 0.0
y0 = 0.0
dydx0 = 0.0

# Define the domain of the solution
a = 0.0
b = L
n = 101   # Number of grid points
h = (b - a) / (n - 1)   # Grid spacing

# Initialize the solution arrays
x = np.linspace(a, b, n)
y = np.zeros(n)

# Define the tridiagonal matrix for the finite difference approximation
A = np.zeros((3, 3))
A[:, 0] = 1.0
A[1:-1, 1] = -(2.0 + (h**2) * (w/(E*I)))
A[1:-1, 2] = 1.0
B = np.zeros(n)
B[1:-1] = -(h**2) * w * x[1:-1]**2 / (2*E*I)

# Apply boundary conditions to the tridiagonal matrix
B[0] = y0
B[-1] = dydx0 * h

print(A)
print(B)

# Solve the tridiagonal system of equations using the Thomas algorithm
C = np.zeros(n-1)
D = np.zeros(n)
C[0] = A[0, 2] / A[0, 1]
D[0] = B[0] / A[0, 1]
for i in range(1, n):
    if i == n-1:
        C[i-1] = 0.0
        D[i] = (B[i] - A[i, 0]*D[i-1]) / A[i, 1]
    else:
        C[i-1] = A[i, 2] / (A[i, 1] - A[i, 0]*C[i-2])
        D[i] = (B[i] - A[i, 0]*D[i-1]) / (A[i, 1] - A[i, 0]*C[i-2])
y[-1] = dydx0*h
y[-2] = D[-1]
for i in range(n-3, -1, -1):
    y[i+1] = C[i]*y[i+2] + D[i]

print(y)

# Plot the solution
plt.plot(x, y*1000, 'r-', lw=2)
plt.xlabel('x (m)')
plt.ylabel('y (mm)')
plt.show()
