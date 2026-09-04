import numpy as np
import matplotlib.pyplot as plt


def trajectory(t, x1_0, x2_0, lam):

    x1 = np.exp(lam * t) * (x1_0 + x2_0 * t)
    x2 = np.exp(lam * t) * x2_0

    return x1, x2


# -------------------------------------------------------
# Initial conditions very close to the origin
# -------------------------------------------------------

initial_conditions = [
    (0.01, 0.005),
    (0.01, 0.01),
    (0.01, 0.02),
    (-0.01, 0.005),
    (-0.01, 0.01),
    (-0.01, 0.02),
    (0.005, -0.005),
    (0.005, -0.01),
    (0.005, -0.02),
]


lambdas = [1, -1]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))


for ax, lam in zip(axes, lambdas):

    # Different time intervals
    if lam == 1:
        t = np.linspace(0, 3, 1500)
    else:
        t = np.linspace(0, 5, 1500)

    # Plot trajectories
    for x1_0, x2_0 in initial_conditions:

        x1, x2 = trajectory(t, x1_0, x2_0, lam)

        ax.plot(x1, x2, linewidth=1.5)

        # Arrow marks
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

    # x2 = 0 eigenvector
    x1_axis = np.linspace(-0.2, 0.2, 300)

    ax.plot(
        x1_axis,
        np.zeros_like(x1_axis),
        linewidth=1.5
    )

    # Origin
    ax.plot(
        0,
        0,
        'ko',
        markersize=5
    )

    # Coordinate axes
    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)

    ax.set_xlabel(r'$x_1$', fontsize=14)
    ax.set_ylabel(r'$x_2$', fontsize=14)

    ax.set_title(
        rf'$\lambda={lam}$',
        fontsize=16
    )

    ax.set_aspect(
        'equal',
        adjustable='box'
    )

    ax.grid(True, alpha=0.2)

    # Limits
    if lam == 1:
        ax.set_xlim(-0.2, 0.2)
        ax.set_ylim(-0.2, 0.2)
    else:
        ax.set_xlim(-0.02, 0.02)
        ax.set_ylim(-0.02, 0.02)


plt.tight_layout()
plt.show()