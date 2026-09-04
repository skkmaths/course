import numpy as np
import matplotlib.pyplot as plt


def trajectory(t, X0, alpha, beta):

    x10, x20 = X0

    x1 = np.exp(alpha * t) * (
        np.cos(beta * t) * x10
        + np.sin(beta * t) * x20
    )

    x2 = np.exp(alpha * t) * (
        -np.sin(beta * t) * x10
        + np.cos(beta * t) * x20
    )

    return x1, x2


# -------------------------------------------------------
# Initial conditions
# -------------------------------------------------------

initial_conditions = [
    (1.0, 0.0),
    (0.0, 1.0),
    (1.0, 1.0),
    (-1.0, 1.0),
    (-1.0, -1.0),
    (1.0, -1.0)
]


# -------------------------------------------------------
# Cases
# -------------------------------------------------------

cases = [
    (0, 1),
    (0, -1),
    (1, 1),
    (1, -1),
    (-1, 1),
    (-1, -1)
]


fig, axes = plt.subplots(2, 3, figsize=(15, 9))


for ax, (alpha, beta) in zip(axes.flat, cases):

    # ---------------------------------------------------
    # Choose time interval
    # ---------------------------------------------------

    if alpha == 0:

        # One complete period
        T = 2 * np.pi / abs(beta)
        t = np.linspace(0, T, 2000)

    else:

        t = np.linspace(0, 4, 2000)

    # ---------------------------------------------------
    # Plot trajectories
    # ---------------------------------------------------

    for X0 in initial_conditions:

        x1, x2 = trajectory(t, X0, alpha, beta)

        ax.plot(x1, x2, linewidth=1.5)

        # ------------------------------------------------
        # Arrow marks
        # ------------------------------------------------

        indices = np.linspace(
            100,
            len(t) - 100,
            5
        ).astype(int)

        for i in indices:

            ax.annotate(
                '',
                xy=(x1[i + 10], x2[i + 10]),
                xytext=(x1[i], x2[i]),
                arrowprops=dict(
                    arrowstyle='->',
                    linewidth=1.2
                )
            )

    # ---------------------------------------------------
    # Equilibrium point
    # ---------------------------------------------------

    ax.plot(0, 0, 'ko', markersize=5)

    # ---------------------------------------------------
    # Coordinate axes
    # ---------------------------------------------------

    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)

    ax.set_xlabel(r'$x_1$', fontsize=13)
    ax.set_ylabel(r'$x_2$', fontsize=13)

    ax.set_title(
        rf'$\alpha={alpha},\quad \beta={beta}$',
        fontsize=14
    )

    ax.set_aspect('equal', adjustable='box')

    # ---------------------------------------------------
    # Limits
    # ---------------------------------------------------

    if alpha == 0:

        ax.set_xlim(-2.2, 2.2)
        ax.set_ylim(-2.2, 2.2)

    elif alpha < 0:

        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)

    else:

        ax.set_xlim(-60, 60)
        ax.set_ylim(-60, 60)


plt.tight_layout()
plt.show()