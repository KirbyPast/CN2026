import math


#   Translation din coeff si v in schema lui horner ca sa putem calcula polinomul intr-un singur punct
def horner(coeffs, v):
    b = coeffs[0]
    for a in coeffs[1:]:
        b = a + b * v
    return b


# Pt fiecare coeficient, returneaza coeficientul * (ordinul - 1), simuland derivarea
def coeficienti_derivata(coeffs):
    n = len(coeffs) - 1
    if n == 0:
        return [0]
    return [coeffs[i] * (n - i) for i in range(n)]

def calc_r(coeffs):
    a0 = abs(coeffs[0])
    A = max(abs(a) for a in coeffs[1:])
    return (a0 + A) / a0

def metoda_newton(coeffs, d_coeffs, x0, epsilon, k_max=1000):
    x = x0
    k = 0
    while True:
        dp = horner(d_coeffs, x)

        if abs(dp) <= epsilon:
            return None, k

        p = horner(coeffs, x)
        dx = p / dp

        x_new = x - dx
        k += 1

        # Conditii de oprire
        if abs(dx) < epsilon:  # Convergenta
            return x_new, k
        if abs(dx) >= 10 ** 8 or k > k_max:  # Divergenta
            return None, k

        x = x_new

def metoda_olver(coeffs, d_coeffs, dd_coeffs, x0, epsilon, k_max=1000):
    x = x0
    k = 0
    while True:
        p = horner(coeffs, x)
        dp = horner(d_coeffs, x)
        ddp = horner(dd_coeffs, x)

        if abs(dp) <= epsilon:
            return None, k

        c_k = (p ** 2 * ddp) / (dp ** 3)

        dx = (p / dp) + 0.5 * c_k

        x_new = x - dx
        k += 1

        if abs(dx) < epsilon:
            return x_new, k
        if abs(dx) >= 10 ** 8 or k > k_max:
            return None, k

        x = x_new

def main():
    P_coeffs = [1.0, -6.0, 11.0, -6.0]
    epsilon = 1e-6

    P_prim = coeficienti_derivata(P_coeffs)
    P_secund = coeficienti_derivata(P_prim)

    # Calculam raza intervalului
    R = calc_r(P_coeffs)
    print(f"Cautam radacinile in intervalul [-{R:.2f}, {R:.2f}]")

    # Generam puncte de start (x0) in intervalul [-R, R]
    # Impartim intervalul in 100 de puncte pentru a gasi cat mai multe radacini
    num_puncte = 100
    puncte_start = [-R + i * (2 * R / num_puncte) for i in range(num_puncte + 1)]

    radacini_newton = []
    radacini_olver = []

    # Memoram rezultatele detaliate pentru fisier
    raport = []

    for x0 in puncte_start:
        # Testam Newton
        r_newton, pasi_n = metoda_newton(P_coeffs, P_prim, x0, epsilon)
        # Testam Olver
        r_olver, pasi_o = metoda_olver(P_coeffs, P_prim, P_secund, x0, epsilon)

        # Filtram radacinile distincte pentru Newton
        if r_newton is not None:
            este_noua = True
            for rn in radacini_newton:
                if abs(r_newton - rn) <= epsilon:
                    este_noua = False
                    break
            if este_noua:
                radacini_newton.append(r_newton)
                raport.append(f"Newton: Radacina {r_newton:.6f} gasita de la x0={x0:.2f} in {pasi_n} pasi.")

        # Filtram radacinile distincte pentru Olver
        if r_olver is not None:
            este_noua = True
            for ro in radacini_olver:
                if abs(r_olver - ro) <= epsilon:
                    este_noua = False
                    break
            if este_noua:
                radacini_olver.append(r_olver)
                raport.append(f"Olver : Radacina {r_olver:.6f} gasita de la x0={x0:.2f} in {pasi_o} pasi.")

    # Afisare pe ecran
    print("\n--- Radacini gasite ---")
    print(f"Metoda Newton (radacini distincte): {[round(r, 4) for r in radacini_newton]}")
    print(f"Metoda Olver (radacini distincte): {[round(r, 4) for r in radacini_olver]}")

    print("\n--- Comparatie Newton vs Olver ---")
    for linie in raport:
        print(linie)

    # Salvare in fisier (doar radacinile distincte se vor memora) [cite: 9, 10]
    with open("radacini_tema7.txt", "w") as f:
        f.write("Radacini distincte - Metoda Newton:\n")
        for r in radacini_newton:
            f.write(f"{r:.6f}\n")

        f.write("\nRadacini distincte - Metoda Olver:\n")
        for r in radacini_olver:
            f.write(f"{r:.6f}\n")

    print("\nRezultatele au fost salvate in 'radacini_tema7.txt'.")


if __name__ == "__main__":
    main()