import numpy as np
import scipy.linalg as la

eps = pow(10, -10)
#if Math.abs(v) > eps: s = 1/v else: nu se poate
n = 5
B = np.random.rand(n, n) #Matrice oarecare pt generare A
A = np.dot(B, B.T)
b = np.random.rand(n)

def rezolvare_lib(A, b):
    P, L, U = la.lu(A)

    x_lib = np.linalg.solve(A, b)
    print("\n--- Descompunerea LU ---")
    print("Matricea P (Permutare):\n", P)
    print("Matricea L:\n", np.round(L, 4))
    print("Matricea U:\n", np.round(U, 4))

    return x_lib

def descompunere_LDLT(A):
    d = np.zeros(n)
    for p in range(n):
        suma_d = 0.0
        for k in range(p):
            # A[p][k] = l_pk de la pasul p
            suma_d += d[k] * (A[p][k] ** 2)

        d[p] = A[p][p] - suma_d

        if abs(d[p]) <= eps:
            print("Eroare, nu se poate face impartirea")
            return None, None

        for i in range(p+1,n):
            suma_l = 0.0
            for k in range(p):
                suma_l += d[k] * A[i][k] * A[p][k]

            A[i][p] = (A[i][p] - suma_l) / d[p]

    return d,A


def calcul_det(d):
    det_A = 1.0
    for valoare in d:
        det_A *= valoare


def rezolva_cholesky(A, d, b, n, eps):
    #Deci, trebuie sa rezolvam pe rand cele 3 ecuatii, incepand cu Lz = b. L e inferior triunghiulara, deci fol subst directa cu formula din pdf
    z = np.zeros(n)
    for i in range(n):
        suma = 0.0
        for j in range(i):
            suma += A[i][j] * z[j]
        z[i] = b[i] - suma

    #Apoi, urmatorul sistem e cu Dy = z
    #Dar D are 0 uri peste tot mai putin diag princ
    y = np.zeros(n)
    for i in range(n):
        if abs(d[i]) <= eps:
            print("Eroare, nu se poate face impartirea")
            return None
        else:
            y[i] = z[i]/d[i]

    x_chol = np.zeros(n)
    for i in range(n-1, -1, -1):
        suma = 0.0
        for j in range(i+1,n):
            suma+=A[j][i] * x_chol[j]
        x_chol[i] = y[i]-suma

    return x_chol

def verificare_norme(A, x_chol, x_lib, b, n):
    v = np.zeros(n)
    for i in range(n):
        suma = 0.0
        for j in range(n):
            #Deci trb sa verificam daca suntem in partea de sus sau partea de jos a matricii.
            #Pentru ca in Aavem acces doar la jumatatea superioara, ca sa aflam elem din jum inferioara trebuie sa dam flip la indici
            if i<=j:
                suma += A[i][j] * x_chol[j]
            else:
                suma += A[j][i] * x_chol[j]
        v[i] = suma

    #Norma 1
    suma_patrate = 0.0
    for i in range(n):
        suma_patrate += (v[i]-b[i]) ** 2
    norma1 = np.sqrt(suma_patrate)

    suma_patrate = 0.0
    for i in range(n):
        suma_patrate += (x_chol[i] - x_lib[i]) ** 2
    norma2 = np.sqrt(suma_patrate)

    return norma1, norma2

if __name__ == "__main__":
    print("-----------------Matricea A--------------")
    print(np.round(A,4))
    x_lib = rezolvare_lib(A, b)
    A_copy = A.copy()
    d, A_copy = descompunere_LDLT(A_copy)
    print("--------------Vectorul d-------------")
    print(np.round(d,4))

    print("--------------Matricea L-------------")
    print(np.round(A_copy,4))

    x_chol = rezolva_cholesky(A_copy, d,b,n,eps)

    print("\n--- Soluția sistemului folosind LDL^T (x_chol) ---")
    print(np.round(x_chol, 12))

    print("\n--- Soluția sistemului de referință (x_lib) ---")
    print(np.round(x_lib, 12))

    norma1, norma2 = verificare_norme(A_copy, x_chol, x_lib, b, n)

    print("\n--- Verificarea Soluției (Norme Euclidiene) ---")
    print(f"||A * x_chol - b||_2 = {norma1:.4e}")
    print(f"||x_chol - x_lib||_2      = {norma2:.4e}")
