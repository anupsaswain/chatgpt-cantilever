import numpy as np
import matplotlib.pyplot as plt

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
n = 101  # Number of grid points
h = (b - a) / (n - 1)   # Grid spacing

# Initialize the solution arrays
x = np.linspace(a, b, n-2)
y = np.zeros(n)
dydx = np.zeros(n)
A = np.zeros((n-2, n-2))
B = np.zeros(n-2)

# Define the tridiagonal matrix for the finite difference approximation
for i in range(0, n-2):
    A[i, i] = 2.0 + (h**2) * (w/(E*I))
    B[i-1] = -(h**2) * y0 * (w/(E*I))

for i in range(1, n-2):
    A[i-1, i] = -1.0
    A[i, i-1] = -1.0


# Apply boundary conditions to the tridiagonal matrix
A[0, 0] = 1.0
B[0] = y0
A[-1, -1] = 1.0
B[-1] = dydx0 * h

# # Solve the tridiagonal system of equations using the Thomas algorithm
# C = np.zeros(n-2)
# D = np.zeros(n-2)
# C[0] = A[0, 1] / A[0, 0]
# D[0] = B[0] / A[0, 0]
# #breakpoint()
# for i in range(1, n-2):
#     D[i] = (B[i] - A[i, i-1]*D[i-1]) / (A[i, i] - A[i, i-1]*C[i-1])
# for i in range(0, n-3):
#     C[i] = A[i, i+1] / (A[i, i] - A[i, i-1]*C[i-1])

# C[n-3] = -1 / (2.0 + (h**2) * (w/(E*I)) + C[n-4])

# y[-1] = dydx0*h
# y[-2] = D[-1]
# for i in range(n-4, -1, -1):
#     y[i+1] = C[i]*y[i+2] + D[i]
# y[0] = y0

y = np.linalg.solve(A, B)

print(A)
print(B)
print(y)
#breakpoint()
# Plot the solution
plt.plot(x, y*1000, 'r-', lw=2, label='FDM')
plt.plot(x, cantilytical(x)*1000, 'b--', lw=4, label='exact')
plt.xlabel('x (m)')
plt.ylabel('y (mm)')
plt.legend()
plt.show()
