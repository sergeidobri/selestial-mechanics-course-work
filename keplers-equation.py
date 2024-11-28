from math import sin, tan, atan, cos, pi, sqrt
import matplotlib.pyplot as plt
import time
from dotenv import load_dotenv
import os


def load_data(path):
    if os.path.exists(path):
        load_dotenv(path)
    else:
        raise FileNotFoundError(f'Файл {path} не был найден.')

    return float(os.getenv('e')), float(os.getenv('frequency'))

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

def main():
    count_time_iterations = 1e5
    path_to_data = 'data.env'
    e, frequency = load_data(path_to_data)
    precision = 1e-5

    frequency /= 24*60*60  # перевод частоты об/сут. -> об/сек.
    period = 1 / frequency
    
    # высчитываем отдельно M и t
    time_start, time_end, time_step = 0, period, period / count_time_iterations

    t_arr = []
    M_arr = []

    while time_start < time_end:
        t_arr.append(time_start)
        M_arr.append(2 * pi * frequency * time_start)
        time_start += time_step

    fig, axs = plt.subplots(2, 2, figsize=(12, 10)) 

    # 1. Метод половинного деления
    print("Начало метода половинного деления.")

    start_time = time.time()

    Nu_arr = []
    E_arr = []
    for index, t in enumerate(t_arr):
        E_arr.append(find_E_half_div(-1e3, 1e3, precision, t, e, frequency, pi))
        Nu_arr.append(computing_Nu(e, E_arr[-1], M_arr[index], pi))

    end_time = time.time()
    execution_time = end_time - start_time
 
    print(f"Время вычисления методом половинного деления: {execution_time} секунд")

    axs[0, 0].plot(t_arr, E_arr, label='E(t)', color='blue', linestyle='-', markersize=0.1)
    axs[0, 0].plot(t_arr, M_arr, label='M(t)', color='red', linestyle='-', markersize=0.1)
    axs[0, 0].plot(t_arr, Nu_arr, label='Nu(t)', color='orange', linestyle='-', markersize=0.1)

    axs[0, 0].set_title('Метод половинного деления')
    axs[0, 0].set_xlabel('Время (t)')
    axs[0, 0].set_ylabel('Значение функции')

    axs[0, 0].grid(True, linestyle='--', alpha=0.5)
    axs[0, 0].legend()
    # 2. Метод золотого сечения
    print("Начало метода золотого сечения.")

    start_time = time.time()

    E_arr = []
    Nu_arr = []

    for index, t in enumerate(t_arr):
        E_arr.append(find_E_golden_ration(-1e3, 1e3, precision, t, e, frequency, pi))
        Nu_arr.append(computing_Nu(e, E_arr[-1], M_arr[index], pi))

    end_time = time.time()
    execution_time = end_time - start_time
 
    print(f"Время вычисления методом золотого сечения: {execution_time} секунд")

    axs[0, 1].plot(t_arr, E_arr, label='E(t)', color='blue', linestyle='-', markersize=0.1)
    axs[0, 1].plot(t_arr, M_arr, label='M(t)', color='red', linestyle='-', markersize=0.1)
    axs[0, 1].plot(t_arr, Nu_arr, label='Nu(t)', color='orange', linestyle='-', markersize=0.1)

    axs[0, 1].set_title('Метод золотого сечения')
    axs[0, 1].set_xlabel('Время (t)')
    axs[0, 1].set_ylabel('Значение функции')

    axs[0, 1].grid(True, linestyle='--', alpha=0.5)
    axs[0, 1].legend()

    # 3. Метод итераций (метод последовательных приближений)
    print("Начало метода итераций (метода последовательных приближений).")

    start_time = time.time()

    E_arr = []
    Nu_arr = []

    for index, t in enumerate(t_arr):
        E_arr.append(find_E_success_approx(M_arr[index], e, precision))
        Nu_arr.append(computing_Nu(e, E_arr[-1], M_arr[index], pi))

    end_time = time.time()
    execution_time = end_time - start_time
 
    print(f"Время вычисления методом итераций (методом последовательных приближений): {execution_time} секунд")

    axs[1, 0].plot(t_arr, E_arr, label='E(t)', color='blue', linestyle='-', markersize=0.1)
    axs[1, 0].plot(t_arr, M_arr, label='M(t)', color='red', linestyle='-', markersize=0.1)
    axs[1, 0].plot(t_arr, Nu_arr, label='Nu(t)', color='orange', linestyle='-', markersize=0.1)

    axs[1, 0].set_title('Метод последовательных приближений')
    axs[1, 0].set_xlabel('Время (t)')
    axs[1, 0].set_ylabel('Значение функции')

    axs[1, 0].grid(True, linestyle='--', alpha=0.5)
    axs[1, 0].legend()

    # 4. Метод Ньютона (метод касательных)

    print("Начало метода Ньютона (метода касательных).")

    start_time = time.time()

    E_arr = []
    Nu_arr = []

    for index, t in enumerate(t_arr):
        E_arr.append(find_E_newton(M_arr[index], e, precision))
        Nu_arr.append(computing_Nu(e, E_arr[-1], M_arr[index], pi))

    end_time = time.time()
    execution_time = end_time - start_time
 
    print(f"Время вычисления методом Ньютона: {execution_time} секунд")

    axs[1, 1].plot(t_arr, E_arr, label='E(t)', color='blue', linestyle='-', markersize=0.1)
    axs[1, 1].plot(t_arr, M_arr, label='M(t)', color='red', linestyle='-', markersize=0.1)
    axs[1, 1].plot(t_arr, Nu_arr, label='Nu(t)', color='orange', linestyle='-', markersize=0.1)

    axs[1, 1].set_title('Метод Ньютона (метод касательных)')
    axs[1, 1].set_xlabel('Время (t)')
    axs[1, 1].set_ylabel('Значение функции')

    axs[1, 1].grid(True, linestyle='--', alpha=0.5)
    axs[1, 1].legend()  

    plt.show()


if __name__ == '__main__':
    main()