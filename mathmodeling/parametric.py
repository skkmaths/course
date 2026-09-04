import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Time interval
t = np.linspace(0, 2*np.pi, 400)

# Parametric curve
x1 =(0.01 + -0.01*t)* np.exp(t)
x2 =-0.01*np.exp(t)

# Create figure with two panels
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

# -------------------------------------------------
# Top panel: t in [0, 2*pi]
# -------------------------------------------------

ax1.set_xlim(0, 2*np.pi)
ax1.set_ylim(-0.5, 0.5)

ax1.set_xlabel(r'$t$')
ax1.set_yticks([])

ax1.set_title(r'Time parameter $t\in[0,2\pi]$')

# Straight line representing the domain
ax1.plot([0, 2*np.pi], [0, 0], 'k-', linewidth=2)

# Moving point on the t-axis
t_point, = ax1.plot([], [], 'o', markersize=10)

# Current t value
t_text = ax1.text(
    0.05, 0.75, '',
    transform=ax1.transAxes
)

# -------------------------------------------------
# Bottom panel: circle in (x1,x2) plane
# -------------------------------------------------

ax2.set_xlim(-1.2, 1.2)
ax2.set_ylim(-1.2, 1.2)

ax2.set_xlabel(r'$x_1$')
ax2.set_ylabel(r'$x_2$')

#ax2.set_title(r'$(x_1,x_2)=(\sin t,\cos t)$')

ax2.set_aspect('equal')
ax2.grid(True)

'''
# Full circle as a reference
theta = np.linspace(0, 2*np.pi, 400)
ax2.plot(np.sin(theta), np.cos(theta),
         '--', alpha=0.3)
'''
# Traced trajectory
circle_line, = ax2.plot([], [], linewidth=2)

# Moving point on circle
circle_point, = ax2.plot([], [], 'o', markersize=10)


# -------------------------------------------------
# Initialization
# -------------------------------------------------

def init():
    t_point.set_data([], [])
    circle_line.set_data([], [])
    circle_point.set_data([], [])
    t_text.set_text('')

    return t_point, circle_line, circle_point, t_text


# -------------------------------------------------
# Animation
# -------------------------------------------------

def update(frame):

    # Current time
    current_t = t[frame]

    # Move point along the straight t-axis
    t_point.set_data([current_t], [0])

    # Update time label
    t_text.set_text(
        r'$t={:.2f}$'.format(current_t)
    )

    # Trace the circle
    circle_line.set_data(
        x1[:frame+1],
        x2[:frame+1]
    )

    # Move point on the circle
    circle_point.set_data(
        [x1[frame]],
        [x2[frame]]
    )

    return t_point, circle_line, circle_point, t_text


# Create animation
ani = FuncAnimation(
    fig,
    update,
    frames=len(t),
    init_func=init,
    interval=20,
    blit=True
)

plt.tight_layout()
plt.show()