from math import sin, tan, atan, cos, pi, sqrt
import os
from dotenv import load_dotenv

def load_info(path_to_file):
    if os.path.exists(path_to_file):
        load_dotenv(path_to_file)
    else:
        raise FileNotFoundError("Файл не был найден")
    
    return float(os.getenv("ECCENTRICITY")), float(os.getenv("PERIOD")), float(os.getenv("GRAVITATION_PARAMETER"))

def computing_Nu(e, E, M_last, pi):
    res = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2))
    if M_last > pi:
        res += 2 * pi
    return res

def find_E_half_div(start, end, precision, t, e, freq, pi):
    E_now = None
    while end - start >= 2 * precision:
        E_now = (start + end) / 2
        f_now = E_now - e * sin(E_now) - 2 * pi * freq * t
        if f_now == 0:
            return E_now
        elif f_now > 0:
            end = E_now
        else:
            start = E_now
    
    if E_now is not None:
        return E_now
    return 'Возникла ошибка при вычислении корня'

def find_E_golden_ration(start, end, precision, t, e, freq, pi):
    E_now = None
    golden_ratio = (5 ** 0.5 + 1) / 2
    while end - start >= 2 * precision:
        E_now = start + (end - start) / golden_ratio
        f_now = E_now - e * sin(E_now) - 2 * pi * freq * t
        if f_now == 0:
            return E_now
        elif f_now > 0:
            end = E_now
        else:
            start = E_now
    
    if E_now is not None:
        return E_now
    return 'Возникла ошибка при вычислении корня'

def find_E_success_approx(M, e, precision):
    E_old = 0
    E_new = M
    while abs(E_new - E_old) >= precision:
        E_old = E_new
        E_new = e * sin(E_old) + M

    return E_new

def find_E_newton(M, e, precision):
    E_old = 0
    E_new = M
    while abs(E_new - E_old) >= precision:
        E_old = E_new
        f_old = E_old - e * sin(E_old) - M
        f_derivative_old = 1 - e * cos(E_old)
        E_new = E_old - (f_old / f_derivative_old)
    
    return E_new
