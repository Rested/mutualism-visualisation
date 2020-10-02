import numpy as np
# from scipy.integrate import odeint
import matplotlib.pyplot as p
import ipywidgets as widgets
from modplot import velovect

initial_r1 = -3
initial_r2 = -10
initial_alpha1 = initial_alpha2 = 0.008
initial_b12 = initial_b21 = 0.4
initial_c1 = initial_c2 = 0.008


def dX_dt(X, t=0, **params):
    animal, plant = X
    return np.array(
        [
            (params["r1"] + (params["b12"] * animal)) * plant
            - (
                    (params["alpha1"] + (params["c1"] * params["b12"] * animal))
                    * (plant ** 2)
            ),
            (params["r2"] + (params["b21"] * plant)) * animal
            - (
                    (params["alpha2"] + (params["c2"] * params["b21"] * plant))
                    * (animal ** 2)
            ),
        ]
    )


def get_fixed_points(**params):
    X_f0 = np.array([0.0, 0.0])
    K1 = params["r1"] / params["alpha1"]
    K2 = params["r2"] / params["alpha2"]
    X_f1 = np.array([K1, 0.0])
    X_f2 = np.array([0.0, K2])
    fixed_points = [X_f0]
    print(dX_dt(list(reversed(X_f1)), **params), dX_dt(list(reversed(X_f2)), **params))
    if all(dX_dt(list(reversed(X_f1)), **params) == 0):
        fixed_points.append(X_f1)
    if all(dX_dt(list(reversed(X_f2)), **params) == 0):
        fixed_points.append(X_f2)

    if params["r1"] > 0 or params["r2"] > 0:
        print(
            f"There are partial extinction points at {fixed_points[1]} and {fixed_points[2]}"
        )
    print("Trivial fixed points:", fixed_points)
    return fixed_points


def get_non_trivial_fixed_points(r1, r2, alpha1, alpha2, b12, b21, c1, c2):
    A = (c2 * b21 * alpha1) + (c1 * b12 * b21)
    B = (alpha1 * alpha2) + (c1 * b12 * r2) - (c2 * b21 * r1) - b12 * b21
    C = -r1 * alpha2 - b12 * r2

    plant_values = np.roots([A, B, C])

    animal_values = list(
        map(lambda x: (r2 + (b21 * x)) / (alpha2 + (c2 * b21 * x)), plant_values)
    )

    non_trivial_fixed_points = list(zip(plant_values, animal_values))

    print("non trivial fixed points", non_trivial_fixed_points)

    for fp in non_trivial_fixed_points:
        print(
            "d/dt close to 0",
            dX_dt(
                reversed(list(fp)),
                r1=r1,
                r2=r2,
                alpha1=alpha1,
                alpha2=alpha2,
                b12=b12,
                b21=b21,
                c1=c1,
                c2=c2,
            ),
        )

    return non_trivial_fixed_points


xmax = ymax = 150
xmin = ymin = 0
p.rcParams["figure.dpi"] = 200


def add_quiver(X1, Y1, normalize_arrows, arrow_style, colour_scheme, **params):
    DX1, DY1 = dX_dt([X1, Y1], **params)  # compute growth rate on the grid
    M = np.hypot(DX1, DY1)  # Norm of the growth rate
    if normalize_arrows:
        M[M == 0] = 1.0  # Avoid zero division errors
        DX1 /= M  # Normalize each arrows
        DY1 /= M
    if arrow_style == "Straight":
        p.quiver(
            Y1,
            X1,
            DX1,
            DY1,
            M,
            pivot="mid",
            cmap=getattr(p.cm, colour_scheme),
            zorder=2,
        )
    elif arrow_style == "Curved":
        velovect(
            p.subplot(),
            Y1.T,
            X1.T,
            DX1.T,
            DY1.T,
            color=M,
            cmap=getattr(p.cm, colour_scheme),
        )
    elif arrow_style == "Stream":
        p.streamplot(
            Y1.T,
            X1.T,
            DX1.T,
            DY1.T,
            color=M,
            cmap=getattr(p.cm, colour_scheme),
            zorder=2,
            maxlength=20,
        )


def add_scatter(fixed_points, non_trivial_fixed_points):
    p.scatter(
        np.array([x[0] for x in fixed_points + non_trivial_fixed_points]),
        np.array([x[1] for x in fixed_points + non_trivial_fixed_points]),
        c="red",
        zorder=4,
    )

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
        add_scatter(triv, non_triv)
    p.xlim(xmin, xmax)
    p.ylim(ymin, ymax)
    p.savefig(f"{image_name}.png", dpi=300)


images = [
    ("heatmap_polarised_intrinsic_growth_rate", {
        "r1": -10,
        "r2": 10,
        "grid_toggle": "Off",
        "num_points": 80,
    }),
    ("heatmap_mutual_dependence", {
        "r1": -10,
        "r2": -10,
        "grid_toggle": "Off",
        "num_points": 80
    }),
    ("heatmap_high_intrinsic_growth", {
        "r1": 10,
        "r2": 10,
        "grid_toggle": "Off",
        "num_points": 80,
    }),
    ("stream_polarised", {
        "r1": 9,
        "r2": -10,
        "arrow_style": "Stream",
        "grid_toggle": "Off",
    }),
    ("stream_polarised_the_other_way_around", {
        "r1": -10,
        "r2": 10,
        "arrow_style": "Stream",
        "grid_toggle": "On",
    }),
    ("stream_mutual_dependence", {
        "r1": -5,
        "r2": -4,
        "arrow_style": "Stream",
        "grid_toggle": "On",
    }),
]


for image_name, params in images:
    save_image(f"landscape_images/{image_name}", **params)