import math
from random import uniform
import time

from Ex3a import my_tan_fraction


c1 = 0.33333333333333333
c2 = 0.133333333333333333
c3 = 0.053968253968254
c4 = 0.0218694885361552


def my_tan_polynomial(x):
    is_negative = False
    if x < 0:
        x = -x
        is_negative = True

    if x > math.pi / 4:
        x = (math.pi / 2) - x
        x_2 = x * x
        x_3 = x_2 * x
        x_5 = x_3 * x * x
        x_7 = x_5 * x * x
        x_9 = x_7 * x * x
        res = 1 / (x + c1 * x_3 + c2 * x_5 + c3 * x_7 + c4 * x_9)
    else:
        x_2 = x * x
        x_3 = x_2 * x
        x_5 = x_3 * x * x
        x_7 = x_5 * x * x
        x_9 = x_7 * x * x
        res = x + c1 * x_3 + c2 * x_5 + c3 * x_7 + c4 * x_9
    return -res if is_negative else res


if __name__ == "__main__":
    diff_f = 0
    diff_max_f = 0
    diff_time_max_f = 0
    diff_max_p = 0
    diff_time_max_p = 0
    precision = pow(10, -5)

    for i in range(10_000):
        x = uniform(-math.pi / 2, math.pi / 2)

        time_before = time.perf_counter()
        my_tan_f = my_tan_fraction(x, precision)
        time_after = time.perf_counter()

        my_tan_time_f = time_after - time_before

        time_before = time.perf_counter()
        my_tan_p = my_tan_polynomial(x)
        time_after = time.perf_counter()

        my_tan_time_p = time_after - time_before

        time_before = time.perf_counter()
        their_tan = math.tan(x)
        time_after = time.perf_counter()

        their_tan_time = time_after - time_before

        diff_time_p = abs(my_tan_time_p - their_tan_time)
        diff_p = abs(my_tan_p - math.tan(x))
        if diff_p > diff_max_p:
            diff_max_p = diff_p
        if diff_time_p > diff_time_max_p:
            diff_time_max_p = diff_time_p

        diff_time_f = abs(my_tan_time_f - their_tan_time)
        diff_f = abs(my_tan_f - math.tan(x))
        if diff_f > diff_max_f:
            diff_max_f = diff_f
        if diff_time_f > diff_time_max_f:
            diff_time_max_f = diff_time_f

    print(
        f"Cea mai mare diferenta (de valoare) dintre functia din librarie si functia noastra (fractii): {diff_max_f}"
    )
    print(
        f"Cea mai mare diferenta (de timp) dintre functia din librarie si functia noastra (fractii): {diff_time_max_f}"
    )
    print(
        f"Cea mai mare diferenta (de valoare) dintre functia din librarie si functia noastra (polinom): {diff_max_p}"
    )
    print(
        f"Cea mai mare diferenta (de timp) dintre functia din librarie si functia noastra (polinom): {diff_time_max_p}"
    )
