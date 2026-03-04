import math
from random import uniform
import time


def my_tan_fraction(x, eps):
    mic = pow(10, -12)
    bj = 0
    f_old = bj

    if f_old == 0:
        f_old = mic

    c_old = f_old
    d_old = 0
    j = 1
    delta = 0

    # start duplicate
    if j == 1:
        aj = x
        bj = 1
    else:
        aj = -pow(x, 2)
        bj = 2 * j - 1

    d_new = bj + aj * d_old

    if d_new == 0:
        d_new = mic

    c_new = bj + (aj / c_old)

    if c_new == 0:
        c_new = mic

    d_new = 1 / d_new
    delta = c_new * d_new
    f_new = delta * f_old
    j = j + 1
    f_old = f_new
    d_old = d_new
    c_old = c_new
    # end duplicate

    while abs(delta - 1) > eps:
        if j == 1:
            aj = x
            bj = 1
        else:
            aj = -pow(x, 2)
            bj = 2 * j - 1

        d_new = bj + aj * d_old

        if d_new == 0:
            d_new = mic

        c_new = bj + (aj / c_old)

        if c_new == 0:
            c_new = mic

        d_new = 1 / d_new
        delta = c_new * d_new
        f_new = delta * f_old
        j = j + 1
        f_old = f_new
        d_old = d_new
        c_old = c_new

    return f_new


if __name__ == "__main__":
    value_diff_max = 0
    time_diff_max = 0
    precision = pow(10, -5)

    for i in range(10_000):
        x = uniform(-math.pi / 2, math.pi / 2)

        time_before = time.perf_counter()
        my_tan = my_tan_fraction(x, precision)
        time_after = time.perf_counter()

        my_time_diff = time_after - time_before

        time_before = time.perf_counter()
        their_tan = math.tan(x)
        time_after = time.perf_counter()

        their_time_diff = time_after - time_before

        value_diff = abs(my_tan - their_tan)
        if value_diff > value_diff_max:
            value_diff_max = value_diff

        time_diff = abs(my_time_diff - their_time_diff)
        if time_diff > time_diff_max:
            time_diff_max = my_time_diff

    print(
        f"Cea mai mare diferenta (de valoare) dintre functia din librarie si functia noastra: {value_diff_max}"
    )
    print(
        f"Cea mai mare diferenta (de timp) dintre functia din librarie si functia noastra: {time_diff_max}"
    )
