import numpy as np


def approximate_gradient(f, x, h=1e-5):
    n = len(x)
    grad = np.zeros(n)

    for i in range(n):
        x_plus_2h = np.copy(x)
        x_plus_2h[i] += 2 * h

        x_plus_h = np.copy(x)
        x_plus_h[i] += h

        x_minus_h = np.copy(x)
        x_minus_h[i] -= h

        x_minus_2h = np.copy(x)
        x_minus_2h[i] -= 2 * h

        f1, f2 = f(x_plus_2h), f(x_plus_h)
        f3, f4 = f(x_minus_h), f(x_minus_2h)
        grad[i] = (-f1 + 8 * f2 - 8 * f3 + f4) / (12 * h)

    return grad


def backtracking_search(f, x, grad):
    eta = 1.0
    beta = 0.8
    p = 1
    grad_norm_sq = np.sum(grad**2)

    while p < 8:
        if f(x - eta * grad) > f(x) - (eta / 2) * grad_norm_sq:
            eta *= beta
            p += 1
        else:
            break

    return eta


def gradient_descent(
    f,
    grad_f,
    x0,
    lr_strategy="const",
    eta_const=1e-3,
    use_approx=False,
    eps=1e-5,
    k_max=30000,
    h=1e-5,
):
    x = np.array(x0, dtype=float)
    k = 0

    while True:
        grad = approximate_gradient(f, x, h) if use_approx else np.array(grad_f(x))
        grad_norm = np.linalg.norm(grad)

        eta = (
            backtracking_search(f, x, grad)
            if lr_strategy == "backtracking"
            else eta_const
        )
        step_norm = eta * grad_norm

        x = x - eta * grad
        k += 1

        if step_norm <= eps:
            return x, k, True
        if k >= k_max or step_norm > 1e10:
            return x, k, False


def sigma(z):
    return 1 / (1 + np.exp(-np.clip(z, -500, 500)))


def f1(x):
    w0, w1 = x[0], x[1]
    return -np.log(1e-15 + 1 - sigma(w0 - w1)) - np.log(1e-15 + sigma(w0 + w1))


def grad_f1(x):
    w0, w1 = x[0], x[1]
    return [sigma(w0 - w1) + sigma(w0 + w1) - 1, sigma(w0 + w1) - sigma(w0 - w1) - 1]


def f2(x):
    return x[0] ** 2 + x[1] ** 2 - 2 * x[0] - 4 * x[1] - 1


def grad_f2(x):
    return [2 * x[0] - 2, 2 * x[1] - 4]


def f3(x):
    return 3 * x[0] ** 2 - 12 * x[0] + 2 * x[1] ** 2 + 16 * x[1] - 10


def grad_f3(x):
    return [6 * x[0] - 12, 4 * x[1] + 16]


def f4(x):
    return x[0] ** 2 - 4 * x[0] * x[1] + 4.5 * x[1] ** 2 - 4 * x[1] + 3


def grad_f4(x):
    return [2 * x[0] - 4 * x[1], -4 * x[0] + 9 * x[1] - 4]


def f5(x):
    return x[0] ** 2 * x[1] - 2 * x[0] * x[1] ** 2 + 3 * x[0] * x[1] + 4


def grad_f5(x):
    return [
        2 * x[0] * x[1] - 2 * x[1] ** 2 + 3 * x[1],
        x[0] ** 2 - 4 * x[0] * x[1] + 3 * x[0],
    ]


if __name__ == "__main__":
    functions = [
        ("F1 (Logistic Loss)", f1, grad_f1, [0.0, 0.0]),
        ("F2", f2, grad_f2, [0.0, 0.0]),
        ("F3", f3, grad_f3, [0.0, 0.0]),
        ("F4", f4, grad_f4, [0.0, 0.0]),
        ("F5", f5, grad_f5, [-0.5, 0.1]),
    ]

    for name, f, grad_f, x0 in functions:
        print(f"--- Testare {name} ---")

        x_an_const, k_an_const, _ = gradient_descent(
            f, grad_f, x0, lr_strategy="const", eta_const=1e-3, use_approx=False
        )
        x_ap_const, k_ap_const, _ = gradient_descent(
            f, grad_f, x0, lr_strategy="const", eta_const=1e-3, use_approx=True
        )

        x_an_bt, k_an_bt, _ = gradient_descent(
            f, grad_f, x0, lr_strategy="backtracking", use_approx=False
        )
        x_ap_bt, k_ap_bt, _ = gradient_descent(
            f, grad_f, x0, lr_strategy="backtracking", use_approx=True
        )

        print(
            f"  Analitic + Const (1e-3):   Iteratii = {k_an_const:<6} Solutie = {np.round(x_an_const, 4)}"
        )
        print(
            f"  Aproximat + Const (1e-3):  Iteratii = {k_ap_const:<6} Solutie = {np.round(x_ap_const, 4)}"
        )
        print(
            f"  Analitic + Backtracking:   Iteratii = {k_an_bt:<6} Solutie = {np.round(x_an_bt, 4)}"
        )
        print(
            f"  Aproximat + Backtracking:  Iteratii = {k_ap_bt:<6} Solutie = {np.round(x_ap_bt, 4)}\n"
        )
