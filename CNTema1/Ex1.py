def precizia_masina():
    m = 0
    u = pow(10, -m)

    while True:
        new_m = m + 1
        new_u = pow(10, -new_m)

        if 1.0 + new_u == 1.0:
            break

        u = new_u
        m = new_m

    return u


if __name__ == "__main__":
    u = precizia_masina()
    print(u)
