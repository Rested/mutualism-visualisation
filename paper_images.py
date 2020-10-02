# from scipy.integrate import odeint
import matplotlib.pyplot as p
import numpy as np

from core.constants import *
from core.maths import get_fixed_points, get_non_trivial_fixed_points
from core.plotting import add_quiver, add_scatter

xmax = ymax = 150
xmin = ymin = 0
p.rcParams["figure.dpi"] = 200


def add_plot(
    plt,
    r1=initial_r1,
    r2=initial_r2,
    alpha1=initial_alpha1,
    alpha2=initial_alpha2,
    c1=initial_c1,
    c2=initial_c2,
    b12=initial_b12,
    b21=initial_b21,
    num_points=20,
    grid_toggle="On",
    fixed_point_toggle="On",
    normalize_arrows="Yes",
    arrow_style="Straight",
    colour_scheme="jet",
):
    x = np.linspace(xmin, xmax, num_points)
    y = np.linspace(ymin, ymax, num_points)

    X1, Y1 = np.meshgrid(x, y)

    # if grid_toggle == "On":
    #     plt.grid()
    triv = get_fixed_points(
        r1=r1, r2=r2, alpha1=alpha1, alpha2=alpha2, b12=b12, b21=b21, c1=c1, c2=c2
    )
    non_triv = get_non_trivial_fixed_points(r1, r2, alpha1, alpha2, b12, b21, c1, c2)
    add_quiver(
        plt,
        X1,
        Y1,
        normalize_arrows == "Yes",
        arrow_style,
        colour_scheme,
        r1=r1,
        r2=r2,
        alpha1=alpha1,
        alpha2=alpha2,
        b12=b12,
        b21=b21,
        c1=c1,
        c2=c2,
    )
    if fixed_point_toggle == "On":
        add_scatter(plt, triv, non_triv)
    plt.set_xlim([xmin, xmax])
    plt.set_ylim([ymin, ymax])
    plt.set_aspect("equal")
    plt.set_title(f"$r_1$ = {r1}, $r_2$ = {r2}")
    plt.set_xlabel("$N_1$")
    plt.set_ylabel("$N_2$")


p.rcParams.update({"font.size": 10})
fig, (ax1, ax2, ax3, ax4) = p.subplots(1, 4, sharey=True, sharex=True, figsize=(20, 6))
examples = [
    {"r1": -10, "r2": -10},
    {"r1": -3, "r2": -10},
    {"r1": 1, "r2": -10},
    {"r1": 1, "r2": 1},
]
for i, plt in enumerate([ax1, ax2, ax3, ax4]):
    add_plot(plt, **{**examples[i], "num_points": 15, "arrow_style": "Curved"})


p.tight_layout()
fig.savefig(f"paper_figure_1.png", dpi=300)
