import numpy as np
import scipy
from scipy import linalg

def jacobi_eigen(A_input, eps=1e-10, kmax=10000):

    A = A_input.astype(float).copy()
    n = A.shape[0]
    U = np.eye(n)  # U = I_n initial (pag. 5)

    for k in range(kmax):
        #Gasim (p,q): indicii elementului maxim nediagonal

        max_val = 0.0
        p, q = 1, 0
        for i in range(1, n):
            for j in range(i):
                if abs(A[i, j]) > max_val:
                    max_val = abs(A[i, j])
                    p, q = i, j

        if max_val < eps:
            break

        # alpha = cot(2*theta) = (a_pp - a_qq) / (2 * a_pq)
        alpha = (A[p, p] - A[q, q]) / (2.0 * A[p, q])

        # t = tg(theta)
        if alpha >= 0:
            t = -alpha + np.sqrt(alpha ** 2 + 1)
        else:
            t = -alpha - np.sqrt(alpha ** 2 + 1)

        c = 1.0 / np.sqrt(1 + t ** 2)
        s = t / np.sqrt(1 + t ** 2)

        # Salvam coloana p si q inainte de modificare
        Ap = A[:, p].copy()
        Aq = A[:, q].copy()

        for j in range(n):
            if j == p or j == q:
                continue
            A[p, j] = c * Ap[j] + s * Aq[j]
            A[j, p] = A[p, j]  # simetrie
            A[q, j] = -s * Ap[j] + c * Aq[j]
            A[j, q] = A[q, j]  # simetrie

        A[p, p] = A[p, p] + t * Ap[q]  # b_pp = a_pp + t*a_pq
        A[q, q] = A[q, q] - t * Ap[q]  # b_qq = a_qq - t*a_pq
        A[p, q] = 0.0
        A[q, p] = 0.0

        Up = U[:, p].copy()
        Uq = U[:, q].copy()
        U[:, p] = c * Up + s * Uq
        U[:, q] = -s * Up + c * Uq

    eigenvalues = np.diag(A)
    return eigenvalues, U, k

def test_jacobi():
    matrices = [
        {
            "name": "Matricea 1 (pag. 9)",
            "A": np.array([[0,0,1],[0,0,1],[1,1,1]], dtype=float),
            "lambda_exact": [-1, 0, 2]
        },
        {
            "name": "Matricea 2 (pag. 10)",
            "A": np.array([[1,1,2],[1,1,2],[2,2,2]], dtype=float),
            "lambda_exact": [0, 2*(1-np.sqrt(2)), 2*(1+np.sqrt(2))]
        },
        {
            "name": "Matricea 3 (pag. 10)",
            "A": np.array([[1,0,1,0],[0,1,0,1],[1,0,1,0],[0,1,0,1]], dtype=float),
            "lambda_exact": [0, 0, 2, 2]
        },
        {
            "name": "Matricea 4 (pag. 10)",
            "A": np.array([[1,2,3,4],[2,3,4,5],[3,4,5,6],[4,5,6,7]], dtype=float),
            "lambda_exact": [0, 0, 2*(4-np.sqrt(21)), 2*(4+np.sqrt(21))]
        },
    ]

    for m in matrices:
        print("=" * 60)
        print(m["name"])
        print("=" * 60)

        A_init = m["A"].copy()
        lambdas, U, iters = jacobi_eigen(A_init)

        print(f"Iteratii: {iters}")
        print(f"Valori proprii Jacobi:  {np.round(lambdas, 6)}")
        print(f"Valori proprii exacte:  {np.round(m['lambda_exact'], 6)}")

        # Verificarea ceruta in tema: ||A_init * U - U * Lambda|| (pag. 1)
        # Folosim lambdas si U exact asa cum le-a returnat Jacobi
        Lambda = np.diag(lambdas)
        diff = A_init @ U - U @ Lambda
        norma = np.linalg.norm(diff)
        print(f"||A_init * U - U * Lambda|| = {norma:.2e}")
        print()


def run_cholesky_series(A, eps=1e-10, kmax=1000):
    print("PARTEA 2: Sirul de matrice prin factorizari Cholesky")

    Ak = A.astype(float).copy()

    for k in range(kmax):
        try:
            L = linalg.cholesky(Ak, lower=True)
        except linalg.LinAlgError:
            print(f"Matricea nu e pozitiv definita la iteratia {k}. Stop.")
            break

        # A(k+1) = L^T * L
        Ak_new = L.T @ L

        # Criteriu oprire: ||A(k+1) - A(k)|| < eps
        if np.linalg.norm(Ak_new - Ak) < eps:
            Ak = Ak_new
            print(f"Converge in {k + 1} iteratii.\n")
            break

        Ak = Ak_new
    else:
        print(f"Nu a convergit in {kmax} iteratii.\n")

    print("Ultima matrice calculata A_final:")
    print(Ak)
    print("\nDiagonala (valorile proprii aproximative):")
    print(np.diag(Ak))


def run_svd(A):
    print("PARTEA 3: SVD (Singular Value Decomposition)")


    p, n = A.shape
    if p < n:
        print("Se face doar pentru p > n!")
        return



    # SVD: A = U * S * V^T  (pag. 8)
    U, sigma, Vt = linalg.svd(A, full_matrices=True)
    V = Vt.T

    print("Valori singulare:", np.round(sigma, 6))

    # Rang = numarul de valori singulare strict pozitive
    eps_rank = 1e-10
    rank = np.sum(sigma > eps_rank)
    print(f"Rangul matricei A: {rank}")

    # Numar de conditionare = sigma_max / sigma_min_pozitiv
    sigma_pos = sigma[sigma > eps_rank]
    cond = sigma_pos.max() / sigma_pos.min()
    print(f"Numar de conditionare k2(A) = {cond}")
    print(f"Numar de conditionare scipy  = {np.linalg.cond(A)}")

    # Pseudoinversa Moore-Penrose: A^I = V * S^I * U^T (pag. 9)
    # Construim S^I ∈ R^{n x p}
    SI = np.zeros((n, p))
    for i in range(rank):
        SI[i, i] = 1.0 / sigma[i]

    AI = V @ SI @ U.T
    print("\nPseudoinversa Moore-Penrose A^I:")
    print(AI)

    # Pseudoinversa least-squares: A^J = (A^T A)^{-1} A^T (pag. 2)
    AJ = np.linalg.inv(A.T @ A) @ A.T

    norma_diff = np.linalg.norm(AI - AJ, ord=1)
    print(f"\n||A^I - A^J||_1 = {norma_diff:.2e}")




if __name__ == "__main__":
    # --- Test jacobi ---
    test_jacobi()


    # --- Test Cholesky series ---
    n = 4
    B = np.random.rand(n, n)
    A_pd = B.T @ B + n * np.eye(n)  # simetrica, pozitiv definita
    print("\n\nMatrice test Cholesky (pozitiv definita, 4x4):")
    print(np.round(A_pd, 3))
    run_cholesky_series(A_pd)

    # --- Test SVD cu p > n ---
    p, n = 5, 3
    A_svd = np.random.randn(p, n)
    print("\n\nMatrice test SVD (5x3):")
    print(np.round(A_svd, 3))
    run_svd(A_svd)

