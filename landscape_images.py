# from scipy.integrate import odeint
import matplotlib.pyplot as p
import numpy as np

from core.constants import *
from core.maths import get_fixed_points, get_non_trivial_fixed_points
from core.plotting import add_quiver, add_scatter

xmax = ymax = 150
xmin = ymin = 0
p.rcParams["figure.dpi"] = 200


def save_image(
    image_name,
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
    p.clf()
    x = np.linspace(xmin, xmax, num_points)
    y = np.linspace(ymin, ymax, num_points)

    X1, Y1 = np.meshgrid(x, y)

    p.title(f"$r_1$ = {r1}, $r_2$ = {r2}")
    p.xlabel("Number of soft landscape elements")
    p.ylabel("Number of hard landscape elements")
    if grid_toggle == "On":
        p.grid(zorder=-10)
    triv = get_fixed_points(
        r1=r1, r2=r2, alpha1=alpha1, alpha2=alpha2, b12=b12, b21=b21, c1=c1, c2=c2
    )
    non_triv = get_non_trivial_fixed_points(r1, r2, alpha1, alpha2, b12, b21, c1, c2)
    add_quiver(
        p,
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
        add_scatter(p, triv, non_triv)
    p.xlim(xmin, xmax)
    p.ylim(ymin, ymax)
    p.savefig(f"{image_name}.png", dpi=300)


images = [
    (
        "heatmap_polarised_intrinsic_growth_rate",
        {"r1": -10, "r2": 10, "grid_toggle": "Off", "num_points": 80,},
    ),
    (
        "heatmap_mutual_dependence",
        {"r1": -10, "r2": -10, "grid_toggle": "Off", "num_points": 80},
    ),
    (
        "heatmap_high_intrinsic_growth",
        {"r1": 10, "r2": 10, "grid_toggle": "Off", "num_points": 80,},
    ),
    (
        "stream_polarised",
        {"r1": 9, "r2": -10, "arrow_style": "Stream", "grid_toggle": "Off",},
    ),
    (
        "stream_polarised_the_other_way_around",
        {"r1": -10, "r2": 10, "arrow_style": "Stream", "grid_toggle": "Off",},
    ),
    (
        "stream_mutual_dependence",
        {"r1": -5, "r2": -4, "arrow_style": "Stream", "grid_toggle": "Off",},
    ),
]

for image_name, params in images:
    save_image(f"landscape_images/{image_name}", **params)
