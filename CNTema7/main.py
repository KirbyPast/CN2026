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

if __name__ == "__main__":
    # Definim lista de polinoame din sectiunea Exemple a PDF-ului
    polinoame = [
        {
            "nume": "P1(x) = (x-1)(x-2)(x-3)",
            "coeffs": [1.0, -6.0, 11.0, -6.0]  # [cite: 82]
        },
        {
            "nume": "P2(x) = (x-2/3)(x-1/7)(x+1)(x-3/2)",
            "coeffs": [42.0, -55.0, -42.0, 49.0, -6.0]  # [cite: 84]
        },
        {
            "nume": "P3(x) = (x-1)(x-1/2)(x-3)(x-1/4)",
            "coeffs": [8.0, -38.0, 49.0, -22.0, 3.0]  # [cite: 87]
        },
        {
            "nume": "P4(x) = (x-1)^2 * (x-2)^2 (Radacini duble)",
            "coeffs": [1.0, -6.0, 13.0, -12.0, 4.0]  # [cite: 90]
        }
    ]

    epsilon = 1e-6 # Precizia pentru oprirea algoritmului (calcul)
    epsilon_distinct = 1e-2 # Precizia mai permisiva pentru filtrarea radacinilor duble (afisare)

    with open("radacini_tema7.txt", "w") as f:
        for idx, pol in enumerate(polinoame, 1):
            P_coeffs = pol["coeffs"]
            print(f"\n{'='*50}\nRezolvare pentru Polinomul {idx}: {pol['nume']}\n{'='*50}")
            f.write(f"\n--- Polinomul {idx}: {pol['nume']} ---\n")

            P_prim = coeficienti_derivata(P_coeffs)
            P_secund = coeficienti_derivata(P_prim)

            R = calc_r(P_coeffs)
            print(f"Cautam radacinile in intervalul [-{R:.2f}, {R:.2f}]")

            num_puncte = 100
            puncte_start = [-R + i * (2 * R / num_puncte) for i in range(num_puncte + 1)]

            radacini_newton = []
            radacini_olver = []
            raport_comparatie = []

            for x0 in puncte_start:
                # Testam metodele
                r_newton, pasi_n = metoda_newton(P_coeffs, P_prim, x0, epsilon)
                r_olver, pasi_o = metoda_olver(P_coeffs, P_prim, P_secund, x0, epsilon)

                # Comparatie pas cu pas pentru radacinile comune gasite de la acelasi x0
                if r_newton is not None and r_olver is not None:
                    # Daca ambele converg spre aceeasi radacina reala
                    if abs(r_newton - r_olver) < epsilon_distinct:
                        raport_comparatie.append(
                            f"x0 = {x0:6.2f} | Rădăcină aprox: {r_newton:8.5f} | Pași Newton: {pasi_n:3} | Pași Olver: {pasi_o:3}"
                        )

                # Filtram radacinile distincte pentru Newton
                if r_newton is not None:
                    este_noua = True
                    for rn in radacini_newton:
                        if abs(r_newton - rn) <= epsilon_distinct:
                            este_noua = False
                            break
                    if este_noua:
                        radacini_newton.append(r_newton)

                # Filtram radacinile distincte pentru Olver
                if r_olver is not None:
                    este_noua = True
                    for ro in radacini_olver:
                        if abs(r_olver - ro) <= epsilon_distinct:
                            este_noua = False
                            break
                    if este_noua:
                        radacini_olver.append(r_olver)

            # Afisare rezultate polinom curent
            print("\nRadacini gasite:")
            print(f"  Newton: {[round(r, 4) for r in radacini_newton]}")
            print(f"  Olver : {[round(r, 4) for r in radacini_olver]}")

            print("\nComparatie nr. de pasi (extras):")
            # Afisam doar primele 10 comparatii ca sa nu umplem ecranul inutil
            for linie in raport_comparatie[:10]:
                print("  " + linie)
            if len(raport_comparatie) > 10:
                print("  ... (restul comparatiilor sunt ascunse)")

            # Salvare in fisier pentru polinomul curent
            f.write("Radacini distincte - Newton:\n")
            for r in radacini_newton:
                f.write(f"{r:.6f}\n")

            f.write("Radacini distincte - Olver:\n")
            for r in radacini_olver:
                f.write(f"{r:.6f}\n")

    print("\nRezultatele complete au fost salvate in 'radacini_tema7.txt'.")
