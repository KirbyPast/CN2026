import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends import BackendFilter, backend_registry


# Nu face neaparat parte din tema
def _should_show_figure() -> bool:
    """True when ``plt.show()`` can display (GUI, Jupyter inline, etc.)."""
    be = matplotlib.get_backend().lower()
    if "inline" in be:
        return True
    interactive = {
        b.casefold() for b in backend_registry.list_builtin(BackendFilter.INTERACTIVE)
    }
    return be in interactive


# Nu face neaparat parte din tema
def _finish_figure(save_path: str) -> None:
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    if _should_show_figure():
        plt.show()
    plt.close()

def horner(a, x):
    d = a[-1]
    for i in range(len(a) - 2, -1, -1):
        d = a[i] + d * x
    return d

def least_squares(x, y, m):
    B = np.zeros((m + 1, m + 1))
    f_vec = np.zeros(m + 1)
    for i in range(m + 1):
        for j in range(m + 1):
            B[i, j] = np.sum(x ** (i + j))
        f_vec[i] = np.sum(y * (x ** i))
    a = np.linalg.solve(B, f_vec)
    return a

def cubic_spline(x, y, da, db):
    n = len(x) - 1
    h = np.diff(x)
    H = np.zeros((n + 1, n + 1))
    rhs = np.zeros(n + 1)

    H[0, 0] = 2 * h[0]
    H[0, 1] = h[0]
    rhs[0] = 6 * ((y[1] - y[0]) / h[0] - da)

    for i in range(1, n):
        H[i, i - 1] = h[i - 1]
        H[i, i] = 2 * (h[i - 1] + h[i])
        H[i, i + 1] = h[i]
        rhs[i] = 6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])

    H[n, n - 1] = h[n - 1]
    H[n, n] = 2 * h[n - 1]
    rhs[n] = 6 * (db - (y[n] - y[n - 1]) / h[n - 1])

    A = np.linalg.solve(H, rhs)
    return A, h

def eval_spline(x_bar, x, y, A, h):
    n = len(x) - 1
    idx = -1
    for i in range(n):
        if x[i] <= x_bar <= x[i + 1]:
            idx = i
            break
    if idx == -1:
        idx = n - 1

    bi = (y[idx + 1] - y[idx]) / h[idx] - h[idx] * (A[idx + 1] - A[idx]) / 6
    ci = (x[idx + 1] * y[idx] - x[idx] * y[idx + 1]) / h[idx] - h[idx] * (x[idx + 1] * A[idx] - x[idx] * A[idx + 1]) / 6

    term1 = ((x_bar - x[idx]) ** 3 * A[idx + 1]) / (6 * h[idx])
    term2 = ((x[idx + 1] - x_bar) ** 3 * A[idx]) / (6 * h[idx])
    
    return term1 + term2 + bi * x_bar + ci

def solve_homework(x0, xn, n, f, da, db, x_bar, m_values, plot_path: str):
    x = np.zeros(n + 1)
    x[0] = x0
    x[n] = xn
    if n > 1:
        x[1:n] = np.sort(np.random.uniform(x0, xn, n - 1))
    y = f(x)
    f_x_bar = f(x_bar)

    print(f"--- Rezultate pentru x_bar = {x_bar} ---")
    print(f"f(x_bar) exact: {f_x_bar:.4f}")

    print("\nMetoda celor mai mici patrate:")
    best_a = None
    best_m = 0
    for m in m_values:
        a = least_squares(x, y, m)
        p_val = horner(a, x_bar)
        err = abs(p_val - f_x_bar)
        sum_err = sum([abs(horner(a, xi) - yi) for xi, yi in zip(x, y)])
        print(f"m={m}: Pm(x_bar)={p_val:.4f}, Eroare={err:.4f}, Suma Erori={sum_err:.4f}")
        if m == m_values[-1]:
            best_a = a
            best_m = m

    print("\nSpline cubic:")
    A, h = cubic_spline(x, y, da, db)
    s_val = eval_spline(x_bar, x, y, A, h)
    print(f"Sf(x_bar)={s_val:.4f}, Eroare={abs(s_val - f_x_bar):.4f}\n")

    x_plot = np.linspace(x0, xn, 100)
    y_exact = f(x_plot)
    y_pm = [horner(best_a, xi) for xi in x_plot]
    y_spline = [eval_spline(xi, x, y, A, h) for xi in x_plot]

    plt.figure(figsize=(10, 6))
    plt.plot(x_plot, y_exact, label='f(x) exact', linewidth=2)
    plt.plot(x_plot, y_pm, label=f'Pm(x) m={best_m}', linestyle='--')
    plt.plot(x_plot, y_spline, label='Sf(x) Spline', linestyle=':')
    plt.scatter(x, y, color='red', label='Noduri', zorder=5)
    plt.axvline(x=x_bar, color='green', linestyle='-', alpha=0.5, label='x_bar')
    plt.legend()
    plt.title('Aproximare interpolare Tema 6')
    plt.grid(True)
    _finish_figure(plot_path)
    print(f"Figura salvata: {plot_path}\n")

if __name__ == '__main__':
    np.random.seed(42)

    def f1(x): return x**4 - 12*x**3 + 30*x**2 + 12
    print("EXEMPLUL 1")
    solve_homework(0, 2, 10, f1, 0, 8, 1.5, [2, 3, 4, 5], "exemplu1.png")

    def f2(x): return x**3 + 3*x**2 - 5*x + 12
    print("EXEMPLUL 2")
    solve_homework(1, 5, 10, f2, 4, 100, 3.0, [2, 3, 4, 5], "exemplu2.png")
