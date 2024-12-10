from keplers_equation import *
import matplotlib.pyplot as plt


def main():
    # Самым эффективным оказался метод последовательных приближений. Используем его.
    e, T, mu = load_info("venus_orbit_elements.env")
    cnt_iter = 100
    precision = 1e-5
    T *= 60 * 60
    freq = 1 / T
    a = (T * T * mu / (4 * pi * pi)) ** (1/3)

    time_start, time_end, time_step = 0, T, T / cnt_iter
    t_arr = []
    Nu_arr = []
    r_arr = []
    velo_r_arr = []
    velo_n_arr = []
    velocity = []
    r_pericentre = 1e10
    r_apocentre = -1

    while time_start < time_end:
        t_arr.append(time_start)

        p_now = a * (1 - e ** 2)
        M_now = 2 * pi * freq * time_start

        Nu_arr.append(computing_Nu(e, find_E_success_approx(M_now, e, precision), M_now, pi))
        r_arr.append(p_now / (1 + e * cos(Nu_arr[-1])))
        velo_r_arr.append((mu / p_now) ** 0.5 * e * sin(Nu_arr[-1]))
        velo_n_arr.append((mu / p_now) ** 0.5 * (1 + e * cos(Nu_arr[-1])))
        velocity.append((velo_r_arr[-1] ** 2 + velo_n_arr[-1] ** 2) ** 0.5)

        if r_arr[-1] > r_apocentre:
            r_apocentre = r_arr[-1]
        if r_arr[-1] < r_pericentre:
            r_pericentre = r_arr[-1]

        time_start += time_step

    fig, axs = plt.subplots(2, 2, figsize=(12, 10)) 
    fig1, fig2, fig3, fig4 = axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]
    
    fig1.plot(t_arr, r_arr, label='r(t)', color='blue', markersize=0.1)
    fig1.set_title('Зависимость радиус-вектора от времени')
    fig1.set_xlabel('Время (t)')
    fig1.set_ylabel('Значение радиус-вектора')
    fig1.grid(True, linestyle='--', alpha=0.5)  
    fig1.legend()

    fig2.plot(t_arr, velo_r_arr, label='v_r(t)', color='red', markersize=0.1)
    fig2.set_title('Зависимость радиальной компоненты скорости от времени')
    fig2.set_xlabel('Время (t)')
    fig2.set_ylabel('Значение радиальной компоненты скорости')
    fig2.grid(True, linestyle='--', alpha=0.5)  
    fig2.legend()

    fig3.plot(t_arr, velo_n_arr, label='v_n(t)', color='green', markersize=0.1)
    fig3.set_title('Зависимость трансверсальной компоненты скорости от времени')
    fig3.set_xlabel('Время (t)')
    fig3.set_ylabel('Значение трансверсальной компоненты скорости')
    fig3.grid(True, linestyle='--', alpha=0.5)  
    fig3.legend()

    fig4.plot(t_arr, velocity, label='v(t)', color='orange', markersize=0.1)
    fig4.set_title('Зависимость скорости от времени')
    fig4.set_xlabel('Время (t)')
    fig4.set_ylabel('Значение скорости')
    fig4.grid(True, linestyle='--', alpha=0.5)  
    fig4.legend()

    plt.show()

if __name__ == '__main__':
    main()
