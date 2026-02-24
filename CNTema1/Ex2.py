import math
from random import uniform
import time


c1=0.33333333333333333
c2=0.133333333333333333
c3=0.053968253968254
c4=0.0218694885361552

def calculate_tan_polynomial(x):
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
        res = 1/(x + c1*x_3 + c2*x_5 + c3*x_7 + c4*x_9)
    else:
        x_2 = x * x
        x_3 = x_2 * x
        x_5 = x_3 * x * x
        x_7 = x_5 * x * x
        x_9 = x_7 * x * x
        res = x + c1*x_3 + c2*x_5 + c3*x_7 + c4*x_9
    return -res if is_negative else res

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
    diff_f = 0
    diff_max_f = -999999
    diff_time_max_f = -99999999
    diff_max_p = -99999
    diff_time_max_p = -99999999
    while i < 10000:
        x = uniform(-math.pi / 2, math.pi / 2)

        time_before = time.perf_counter()
        my_tan_f = calculate_tan_fraction(x, pow(10, -5))
        time_after = time.perf_counter()

        my_tan_time_f = time_after - time_before

        time_before = time.perf_counter()
        my_tan_p = calculate_tan_polynomial(x)
        time_after = time.perf_counter()

        my_tan_time_p = time_after - time_before

        time_before = time.perf_counter()
        their_tan = math.tan(x)
        time_after = time.perf_counter()

        their_tan_time = time_after - time_before
        diff_time_f = abs(my_tan_time_f - their_tan_time)
        diff_f = abs(my_tan_f - math.tan(x))
        diff_time_p = abs(my_tan_time_p - their_tan_time)
        diff_p = abs(my_tan_p - math.tan(x))
        if diff_p > diff_max_p:
            diff_max_p = diff_p
        if diff_time_p > diff_time_max_p:
            diff_time_max_p = diff_time_p

        if diff_f > diff_max_f:
            diff_max_f = diff_f
        if diff_time_f > diff_time_max_f:
            diff_time_max_f = diff_time_f
        i = i + 1
    print("Cea mai mare diferenta dintre functia din librarie si functia noastra (fractii):" + diff_max_f.__str__())
    print("Cea mai mare diferenta (de timp) dintre functia din librarie si functia noastra (fractii):" + diff_time_max_f.__str__())
    print("Cea mai mare diferenta dintre functia din librarie si functia noastra (polinom):" + diff_max_p.__str__())
    print("Cea mai mare diferenta (de timp) dintre functia din librarie si functia noastra(polinom):" + diff_time_max_p.__str__())