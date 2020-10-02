from core.maths import dX_dt
import numpy as np
from modplot import velovect
from matplotlib.cm import get_cmap


def add_quiver(p, X1, Y1, normalize_arrows, arrow_style, colour_scheme, **params):
    DX1, DY1 = dX_dt([X1, Y1], **params)  # compute growth rate on the grid
    M = np.hypot(DX1, DY1)  # Norm of the growth rate
    if normalize_arrows:
        M[M == 0] = 1.0  # Avoid zero division errors
        DX1 /= M  # Normalize each arrows
        DY1 /= M
    if arrow_style == "Straight":
        p.quiver(
            Y1, X1, DX1, DY1, M, pivot="mid", cmap=get_cmap(colour_scheme), zorder=2,
        )
    elif arrow_style == "Curved":
        try:
            sp = p.subplot()
        except AttributeError:
            sp = p
        velovect(
            sp,
            Y1.T,
            X1.T,
            DX1.T,
            DY1.T,
            color=M,
            cmap=get_cmap(colour_scheme),
            # arrowsize=20
        )
    elif arrow_style == "Stream":
        p.streamplot(
            Y1.T,
            X1.T,
            DX1.T,
            DY1.T,
            color=M,
            cmap=get_cmap(colour_scheme),
            zorder=2,
            maxlength=20,
        )


def add_scatter(p, fixed_points, non_trivial_fixed_points):
    p.scatter(
        np.array([x[0] for x in fixed_points + non_trivial_fixed_points]),
        np.array([x[1] for x in fixed_points + non_trivial_fixed_points]),
        c="red",
        zorder=4,
    )
