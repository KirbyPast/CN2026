import math
from random import uniform
import time

def calculate_tan_fraction(x,eps):
    mic = pow(10, -12)
    bj = 0
    f_old = bj
    if f_old == 0:
        f_old = mic
    c_old = f_old
    d_old = 0
    j = 1
    delta = 0
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

if __name__=="__main__":
    i = 1
    diff = 0
    diff_max = -999999
    diff_time_max = -99999999
    while i < 10000:
        x = uniform(-math.pi / 2, math.pi / 2)

        time_before = time.perf_counter()
        my_tan = calculate_tan_fraction(x,pow(10, -5))
        time_after = time.perf_counter()

        my_tan_time = time_after-time_before

        time_before = time.perf_counter()
        their_tan = math.tan(x)
        time_after = time.perf_counter()

        their_tan_time = time_after - time_before
        diff_time = abs(my_tan_time - their_tan_time)
        diff = abs(my_tan - math.tan(x))
        if diff > diff_max:
            diff_max = diff
        if diff_time > diff_time_max:
            diff_time_max = diff_time
        i = i + 1
    print("Cea mai mare diferenta dintre functia din librarie si functia noastra:" + diff_max.__str__())
    print("Cea mai mare diferenta (de timp) dintre functia din librarie si functia noastra:" + diff_time_max.__str__())
