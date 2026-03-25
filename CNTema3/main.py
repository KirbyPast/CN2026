import math

import numpy as np

eps = pow(10,-9)


def init_data(n):
    A = np.random.rand(n, n) * 15
    s = np.random.rand(n) * 15
    return A, s

def calc_b(A,s):
    n=s.size
    b = np.zeros(n)
    for i in range(n):
        for j in range(n):
            b[i] += s[j] * A[i][j]

    return b

def desc_QR_householder(A_start, b_start):
    n = len(A_start)

    A = np.copy(A_start)
    b = np.copy(b_start)
    Q_barat = np.eye(n)

    for r in range(n-1):
        sigma = 0
        for j in range(r,n):
            sigma += A[j][r]*A[j][r]

        if sigma <= eps:
            continue

        k = math.sqrt(sigma)
        if A[r][r] > 0:
            k = -k

        beta = sigma - k*A[r][r]

        u = np.zeros(n)
        u[r] = A[r][r] - k
        for i in  range(r+1, n):
            u[i] = A[i][r]
        # Transformam A
        for j in range(r+1, n):
            gamma = 0.0
            for i in range(r, n):
                gamma += u[i] * A[i][j]
            gamma /= beta
            for i in range(r,n):
                A[i][j] = A[i][j] - gamma * u[i]

        A[r][r] = k
        for i in range(r+1,n):
            A[i][r] = 0
        #Transformam b
        gamma = 0.0
        for i in range(r,n):
            gamma += u[i]*b[i]
        gamma /= beta

        for i in range(r,n):
            b[i] = b[i] - gamma * u[i]
        #Transformam Q barat

        for j in range(n):
            gamma = 0.0
            for i in range(r,n):
                gamma += u[i] * Q_barat[i][j]
            gamma /= beta
            for i in range(r,n):
                Q_barat[i][j] -= gamma * u[i]

    return A,b,Q_barat


def substitutie_inversa(A, b):
    n = len(b)
    x = np.zeros(n)

    for i in range(n - 1, -1, -1):
       if abs(A[i][i]) <= eps:
           print(f"Eroare: Matricea este singulara (elementul A[{i}][{i}] este ~0).")
           return None

       suma = 0.0
       for j in range(i + 1, n):
          suma += A[i][j] * x[j]

       x[i] = (b[i] - suma) / A[i][i]

    return x

def calc_norma(x1, x2):
    n = len(x1)
    suma_patrate = 0.0

    for i in range(n):
        diferenta = x1[i] - x2[i]
        suma_patrate += diferenta * diferenta
    return math.sqrt(suma_patrate)

def calc_norma_matrice(M1, M2):
    n = len(M1)
    suma_patrate = 0.0
    for i in range(n):
        for j in range(n):
            diferenta = M1[i][j] - M2[i][j]
            suma_patrate += diferenta * diferenta
    return math.sqrt(suma_patrate)

def inversare_QR(R, Q):
    n = len(R)
    A_inv = np.zeros((n,n))

    for i in range(n):
        if abs(R[i][i]) < eps:
            print("Eroare, inversa nu poate fi calculata!")
            return None

    for j in range(n):
        b_j = np.zeros(n)
        for i in range(n):
            b_j[i] = Q[i][j]

        x_col = substitutie_inversa(R,b_j)

        for i in range(n):
            A_inv[i][j] = x_col[i]

    return A_inv


if __name__ == "__main__":
    np.set_printoptions(precision=7)
    # #6 din tema
    n = int(input("Introduceti marimea datelor: "))
    A, s = init_data(n)
    print("Matricea A:")
    print(A)
    print("Vectorul s:")
    print(s)
    # #1 din tema
    b = calc_b(A, s)
    print("Vectorul b:")
    print(b)
    # #2 din tema
    A_nou, b_nou, Q_t = desc_QR_householder(A,b)
    print("Matricea A (R): (noua)")
    print(A_nou)
    print("Vectorul b: (nou)")
    print(b_nou)
    print("Matricea Q: ")
    print(Q_t.T)
    print("Q * Qt")
    print(np.dot(Q_t,Q_t.T))
    # #3 din tema
    x_householder = substitutie_inversa(A_nou, b_nou)

    Q_lib, R_lib = np.linalg.qr(A)
    b_lib = np.dot(Q_lib.T, b)
    x_QR = substitutie_inversa(R_lib, b_lib)

    print("x-ul nostru:")
    print(x_householder)
    print("x-ul librariei:")
    print(x_QR)
    if x_householder is not None:
        print("Eroare intre x house si x qr: ")
        print(calc_norma(x_QR, x_householder))

        # #4 din tema
        print("Eroare intre A_init*x_householder si b init: ")
        print(calc_norma(np.dot(A, x_householder), b)) #prima marja, A init * x - b init
        print("Eroare intre A_init*x_qr si b init: ")
        print(calc_norma(np.dot(A, x_QR), b)) #A doua, A init * xqr - b init
        print("Eroare intre x_householder si s: ")
        print(calc_norma(x_householder, s) / calc_norma(s, np.zeros(n)))
        print("Eroare intre x_qr si s: ")
        print(calc_norma(x_QR, s) / calc_norma(s, np.zeros(n)))

    # #5 din tema
    A_inv_householder = inversare_QR(A_nou, Q_t)
    A_inv_lib = np.linalg.inv(A)
    if A_inv_householder is not None:
        print("Eroare intre A_inv_house si A_inv_lib: ")
        print(calc_norma_matrice(A_inv_lib, A_inv_householder))

