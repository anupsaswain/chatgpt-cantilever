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
n = 101   # Number of grid points
h = (b - a) / (n - 1)   # Grid spacing

# Initialize the solution arrays
x = np.linspace(a, b, n)
y = np.zeros(n)
dydx = np.zeros(n)

# Define the function for the right-hand side of the differential equation
def f(x, y, dydx):
    return np.array([dydx, -w/(E*I)])

# Apply fourth-order Runge-Kutta method to solve the differential equation
for i in range(1, n):
    k1 = h * f(x[i-1], y[i-1], dydx[i-1])
    k2 = h * f(x[i-1] + h/2, y[i-1] + k1[0]/2, dydx[i-1] + k1[1]/2)
    k3 = h * f(x[i-1] + h/2, y[i-1] + k2[0]/2, dydx[i-1] + k2[1]/2)
    k4 = h * f(x[i-1] + h, y[i-1] + k3[0], dydx[i-1] + k3[1])
    y[i] = y[i-1] + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6
    dydx[i] = dydx[i-1] + (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6

# Apply boundary conditions
y[0] = y0
dydx[0] = dydx0

# Plot the solution
plt.plot(x, y*1000, 'r-', lw=2, label='FDM')
plt.plot(x, cantilytical(x)*1000, 'b--', lw=4, label='exact')
plt.xlabel('x (m)')
plt.ylabel('y (mm)')
plt.legend()
plt.show()
