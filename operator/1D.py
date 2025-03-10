# %%
import numpy as np
import matplotlib.pyplot as plt

def periodic_diff(diff, L):
    """
    Map a distance into the interval [-L/2, L/2] for periodic boundaries.
    """
    return (diff + L/2) % L - L/2

def peridynamic_first_derivative(u, dx, delta, L):
    """
    Approximate the first derivative u_x using a peridynamic operator.
    
    Parameters:
      u     : array of solution values at grid points
      dx    : grid spacing
      delta : horizon for the nonlocal interaction
      L     : length of the periodic domain
      
    Returns:
      d1 : array approximating u_x
    """
    N = len(u)
    d1 = np.zeros_like(u)
    # Determine number of points within the horizon
    horizon_points = int(np.ceil(delta/dx))
    
    for i in range(N):
        sum_val = 0.0
        x_i = i * dx
        # Sum contributions from neighbors within the horizon
        for j in range(i - horizon_points, i + horizon_points + 1):
            j_mod = j % N  # periodic index
            x_j = j_mod * dx
            xi = periodic_diff(x_j - x_i, L)
            if np.abs(xi) <= delta and j_mod != i:
                sum_val += xi * (u[j_mod] - u[i])
        # Normalize the sum to approximate the derivative
        d1[i] = (3.0 / (2 * delta**3)) * sum_val * dx
    return d1

# %%
# ---------------------------
# Problem Setup
# ---------------------------
# Domain parameters
L = 1.0         # Length of the domain
N = 200         # Number of grid points
dx = L / N
x = np.linspace(0, L - dx, N)

# Horizon for the peridynamic operator (e.g., 5 grid spacings)
delta = 10 * dx

# Advection parameter
c = 1.0         # Advection speed

# Time-stepping parameters
T_final = 0.1
dt = 0.00010
nsteps = int(T_final / dt)

# Initial condition: Gaussian pulse centered at x=0.5
u0 = np.exp(-200 * (x - 0.5)**2)
u = u0.copy()

# ---------------------------
# Time Integration (Explicit Euler)
# ---------------------------
for step in range(nsteps):
    d1 = peridynamic_first_derivative(u, dx, delta, L)
    # Update using u_t + c*u_x = 0  => u_t = - c*u_x
    u = u - dt * c * d1

# %%
# ---------------------------
# Plotting the Results
# ---------------------------
plt.figure(figsize=(8, 4))
plt.plot(x, u0, label='Initial')
plt.plot(x, u, label='Final')
plt.xlabel('x')
plt.ylabel('u')
plt.legend()
plt.title('1D Advection Equation via Peridynamic Differential Operator')
plt.show()

# %%
