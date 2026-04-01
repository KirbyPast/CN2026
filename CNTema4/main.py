def read_vector(filepath):
    with open(filepath, "r") as f:
        return [float(line.strip()) for line in f if line.strip()]


def solve_sparse_system(file_prefix="1", epsilon=1e-6, k_max=10000):
    try:
        d0 = read_vector(f"date/d0_{file_prefix}.txt")
        d1 = read_vector(f"date/d1_{file_prefix}.txt")
        d2 = read_vector(f"date/d2_{file_prefix}.txt")
        b = read_vector(f"date/b_{file_prefix}.txt")
    except FileNotFoundError as e:
        print(f"Eroare la citirea fisierelor: {e}")
        return

    # 1. Dimensiunea sistemului
    n = len(d0)
    if len(b) != n:
        print("Eroare: Dimensiunile vectorilor d0 si b nu corespund.")
        return

    # 2. Determinarea relatiilor pentru p si q
    p = n - len(d1)
    q = n - len(d2)

    print(f"--- Sistemul {file_prefix} ---")
    print(f"Dimensiunea sistemului (n): {n}")
    print(f"Diagonala p: {p}")
    print(f"Diagonala q: {q}")

    # 3. Verificarea elementelor nenule pe diagonala principala
    for val in d0:
        if abs(val) <= epsilon:
            print(
                "Eroare: Exista elemente nule (sau < epsilon) pe diagonala principala."
            )
            return
    print("Nu exista elemente nule pe diagonala principala.")

    # 4. Metoda Gauss-Seidel
    x = [0.0] * n
    converged = False

    for k in range(k_max):
        x_prev = x[:]

        for i in range(n):
            sum_val = b[i]

            # Diagonalele p
            if i >= p:
                #(j < i)
                sum_val -= d1[i - p] * x[i - p]
            if i + p < n:
                #(j > i)
                sum_val -= d1[i] * x_prev[i + p]

            # Diagonalele q
            if i >= q:
                #(j < i)
                sum_val -= d2[i - q] * x[i - q]
            if i + q < n:
                #(j > i)
                sum_val -= d2[i] * x_prev[i + q]

            x[i] = sum_val / d0[i]
        # Criteriul de oprire
        dx = max(abs(x[i] - x_prev[i]) for i in range(n))

        if dx < epsilon:
            converged = True
            print(f"Convergenta atinsa in {k + 1} iteratii.")
            break
        elif dx > 1e10:
            print("Sistemul diverge.")
            return

    if not converged:
        print("Numarul maxim de iteratii a fost atins fara convergenta.")
        return

    # 5. Calculul y = A * x_GS
    y = [0.0] * n
    for i in range(n):
        y[i] = d0[i] * x[i]
        if i >= p:
            y[i] += d1[i - p] * x[i - p]
        if i + p < n:
            y[i] += d1[i] * x[i + p]
        if i >= q:
            y[i] += d2[i - q] * x[i - q]
        if i + q < n:
            y[i] += d2[i] * x[i + q]

    # 6. Calculul normei ||A*x_GS - b||_inf
    err_norm = max(abs(y[i] - b[i]) for i in range(n))
    print(f"Norma ||A*x_GS - b||_inf: {err_norm}")

    # Extra: comparatie cu solutia data
    x_exact = get_exact_solution(file_prefix, n)
    max_error = max(abs(x[i] - x_exact[i]) for i in range(n))
    print(f"Eroarea maxima fata de solutia exacta: {max_error}\n")


def get_exact_solution(file_prefix, n):
    idx = int(file_prefix)
    if idx == 1:
        return [2.0 / 3.0 for _ in range(n)]
    elif idx == 2:
        return [i * (4.0 / 7.0) for i in range(n)]
    elif idx == 3:
        return [2.0 for _ in range(n)]
    elif idx == 4:
        return [float(i) for i in range(n)]
    elif idx == 5:
        return [i * 0.4 for i in range(n)]
    return [0.0] * n


if __name__ == "__main__":
    for i in range(5):
        solve_sparse_system(file_prefix=str(i + 1), epsilon=1e-7)
