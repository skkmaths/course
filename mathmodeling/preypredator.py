import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------
# Parameters
# -------------------------------------------------
a = 1
b = 1
m = 1
n = 1

# -------------------------------------------------
# Grid
# -------------------------------------------------
x = np.linspace(0, 4, 25)
y = np.linspace(0, 4, 25)

X, Y = np.meshgrid(x, y)

# -------------------------------------------------
# Vector field
# dx/dt = (a - by)x
# dy/dt = (-m + nx)y
# -------------------------------------------------
dX = (a - b * Y) * X
dY = (-m + n * X) * Y

# Normalize arrows
magnitude = np.sqrt(dX**2 + dY**2)

U = dX / (magnitude + 1e-12)
V = dY / (magnitude + 1e-12)

# -------------------------------------------------
# Plot direction field
# -------------------------------------------------
plt.figure(figsize=(8, 7))

plt.quiver(
    X, Y, U, V,
    angles='xy',
    scale_units='xy',
    scale=12,
    width=0.0025
)

# -------------------------------------------------
# Nullclines
# dx/dt = 0: x = 0 or y = a/b
# dy/dt = 0: y = 0 or x = m/n
# -------------------------------------------------

plt.axhline(
    a / b,
    linestyle='--',
    linewidth=1.5,
    label=r'$\dot{x}=0:\quad y=a/b$'
)

plt.axvline(
    m / n,
    linestyle='--',
    linewidth=1.5,
    label=r'$\dot{y}=0:\quad x=m/n$'
)

# Interior equilibrium
xe = m / n
ye = a / b

plt.plot(
    xe, ye,
    'o',
    markersize=7,
    label=r'Equilibrium $(m/n,a/b)$'
)

# -------------------------------------------------
# Labels and formatting
# -------------------------------------------------
plt.xlabel(r'$x$', fontsize=15)
plt.ylabel(r'$y$', fontsize=15)

plt.xlim(0, 4)
plt.ylim(0, 4)

plt.grid(alpha=0.3)
plt.legend(fontsize=11)

plt.tight_layout()
plt.show()