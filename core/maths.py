import numpy as np


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
