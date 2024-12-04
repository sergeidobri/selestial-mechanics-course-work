from keplers_equation import *
import matplotlib.pyplot as plt
import time


def main():
    count_time_iterations = 1e5  # количество итераций
    e, period, mu = load_info("selestial_mechanics/venus_orbit_elements.env")
    period *= 3600
    precision = 1e-10

    # frequency /= 24*60*60  # перевод частоты об/сут. -> об/сек.
    frequency = 1 / period  # частота обращения
    
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

    fig1 = axs[0, 0]

    fig1.plot(t_arr, E_arr, label='E(t)', color='blue', markersize=0.1)
    fig1.plot(t_arr, M_arr, label='M(t)', color='red', markersize=0.1)
    fig1.plot(t_arr, Nu_arr, label='Nu(t)', color='orange', markersize=0.1)

    fig1.set_title('Метод половинного деления')
    fig1.set_xlabel('Время (t)')
    fig1.set_ylabel('Значение функции')

    fig1.grid(True, linestyle='--', alpha=0.5)
    fig1.legend()

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

    fig2 = axs[0, 1]

    fig2.plot(t_arr, E_arr, label='E(t)', color='blue', markersize=0.1)
    fig2.plot(t_arr, M_arr, label='M(t)', color='red', markersize=0.1)
    fig2.plot(t_arr, Nu_arr, label='Nu(t)', color='orange', markersize=0.1)

    fig2.set_title('Метод золотого сечения')
    fig2.set_xlabel('Время (t)')
    fig2.set_ylabel('Значение функции')

    fig2.grid(True, linestyle='--', alpha=0.5)
    fig2.legend()

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

    fig3 = axs[1, 0]

    fig3.plot(t_arr, E_arr, label='E(t)', color='blue', markersize=0.1)
    fig3.plot(t_arr, M_arr, label='M(t)', color='red', markersize=0.1)
    fig3.plot(t_arr, Nu_arr, label='Nu(t)', color='orange', markersize=0.1)

    fig3.set_title('Метод последовательных приближений')
    fig3.set_xlabel('Время (t)')
    fig3.set_ylabel('Значение функции')

    fig3.grid(True, linestyle='--', alpha=0.5)
    fig3.legend()

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

    fig4 = axs[1, 1]

    fig4.plot(t_arr, E_arr, label='E(t)', color='blue', markersize=0.1)
    fig4.plot(t_arr, M_arr, label='M(t)', color='red', markersize=0.1)
    fig4.plot(t_arr, Nu_arr, label='Nu(t)', color='orange', markersize=0.1)

    fig4.set_title('Метод Ньютона (метод касательных)')
    fig4.set_xlabel('Время (t)')
    fig4.set_ylabel('Значение функции')

    fig4.grid(True, linestyle='--', alpha=0.5)
    fig4.legend()  

    plt.show()

if __name__ == '__main__':
    main()