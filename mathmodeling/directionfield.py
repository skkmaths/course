import numpy as np
import matplotlib.pyplot as plt

# ODE: y' = sin(t)
def f(t, y):
    #return np.sin(t)
    return y-t*t

# Initial condition
t0 = 0
y0 = 1 # -3/2

# --------------------------------------------------
# 1. Create the direction field
# --------------------------------------------------

t_values = np.linspace(-2*np.pi, 2*np.pi, 25)
y_values = np.linspace(-3, 3, 25)

T, Y = np.meshgrid(t_values, y_values)

# Slopes
S = f(T, Y)

# Normalize arrows so that they have comparable lengths
U = np.ones_like(S)
V = S

N = np.sqrt(U**2 + V**2)
U = U / N
V = V / N

plt.figure(figsize=(10, 6))

# Direction field
plt.quiver(T, Y, U, V, angles='xy', pivot='mid',
           color='gray', alpha=0.7)

# --------------------------------------------------
# 2. Follow the solution forward from t = 0
# --------------------------------------------------

h = 0.01
t_forward = np.arange(t0, 2*np.pi + h, h)

y_forward = np.zeros_like(t_forward)
y_forward[0] = y0

for i in range(len(t_forward)-1):
    y_forward[i+1] = y_forward[i] + h*f(t_forward[i], y_forward[i])

# --------------------------------------------------
# 3. Follow the solution backward from t = 0
# --------------------------------------------------

t_backward = np.arange(t0, -2*np.pi - h, -h)

y_backward = np.zeros_like(t_backward)
y_backward[0] = y0

for i in range(len(t_backward)-1):
    y_backward[i+1] = y_backward[i] - h*f(t_backward[i], y_backward[i])

# Reverse backward solution for plotting
t_backward = t_backward[::-1]
y_backward = y_backward[::-1]

# Combine forward and backward solutions
t_solution = np.concatenate((t_backward, t_forward[1:]))
y_solution = np.concatenate((y_backward, y_forward[1:]))

# --------------------------------------------------
# 4. Plot the solution
# --------------------------------------------------

plt.plot(t_solution, y_solution, 'r', linewidth=2.5,
         label='Numerical solution')

plt.plot(t0, y0, 'ko', markersize=7,
         label=r'Initial condition $(0,-3/2)$')

plt.xlabel('$t$')
plt.ylabel('$y$')
plt.title(r"Direction field for $y'=\sin(t)$")
plt.xlim(-2*np.pi, 2*np.pi)
plt.ylim(-3, 3)
plt.grid(True)
plt.legend()

plt.savefig("sol.pdf", bbox_inches="tight")
plt.show()