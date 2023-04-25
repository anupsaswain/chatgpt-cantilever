import numpy as np
import matplotlib.pyplot as plt

# Define the properties of the cantilever
L = 1.0     # Length of cantilever
E = 1.0e7   # Young's modulus
I = 1.0e-4  # Moment of inertia
w = 100.0   # Tip load

# A's added definition of analytical solution for comparison
def cantilytical(x):
    return (w*(x**3)/(6*E*I) - w*L*x*x/(2*E*I))

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
dydx = np.zeros(n)

# Apply finite difference approximation to solve the differential equation
for i in range(1, n):
    y[i] = y[i-1] + h*dydx[i-1]
    dydx[i] = dydx[i-1] + h*(-w/(E*I))
    
# Apply boundary conditions
y[0] = y0
dydx[0] = dydx0

print(y)

# Plot the solution
plt.plot(x, y*1000, 'r-', lw=2, label='FDM')
plt.plot(x, cantilytical(x)*1000, 'b--', lw=4, label='exact')
plt.xlabel('x (m)')
plt.ylabel('y (mm)')
plt.legend()
plt.show()

# yay!! :D