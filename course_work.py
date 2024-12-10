import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from keplers_equation import find_E_newton, computing_Nu
from dotenv import load_dotenv
import os.path

# информация об орбите взята отсюда:
# https://ntrs.nasa.gov/api/citations/19930020400/downloads/19930020400.pdf

def load_orbit_info(path_to_config):
    if os.path.exists(path_to_config):
        load_dotenv(path_to_config)
    else:
        raise FileNotFoundError("Файл не был найден")
    
    return (
        float(os.getenv("ECCENTRICITY")), 
        float(os.getenv("PERIOD")), 
        float(os.getenv("GRAVITATION_PARAMETER")),
        float(os.getenv("ARGUMENT_OF_PERIGEE")),
        float(os.getenv("LONGITUDE_OF_ASCENDING_NODE")),
        float(os.getenv("INCLINATION"))
        )


def deg2rad(degrees):
    return degrees * np.pi / 180

def orbital_elements_to_cartesian(a, e, i, Omega, omega, T, num_points=500):
    t = np.linspace(0, T, num_points)

    nu = np.zeros_like(t)

    for index, time_now in enumerate(t):
        M = 2 * np.pi * time_now / T
        nu[index] = computing_Nu(e, find_E_newton(M, e, 1e-5), M, np.pi)

    r = (a * (1 - e**2)) / (1 + e * np.cos(nu))
    
    x_orb = r * np.cos(nu)
    y_orb = r * np.sin(nu)
    z_orb = np.zeros_like(x_orb)  # массив из num_points нулей
    
    i_rad = deg2rad(i)              # наклонение
    Omega_rad = deg2rad(Omega)      # долгота восходящего узла
    omega_rad = deg2rad(omega)      # аргумент перигея
    
    R_z_Omega = np.array([
        [np.cos(Omega_rad), -np.sin(Omega_rad), 0],
        [np.sin(Omega_rad), np.cos(Omega_rad), 0],
        [0, 0, 1]
    ])
    
    R_x_i = np.array([
        [1, 0, 0],
        [0, np.cos(i_rad), -np.sin(i_rad)],
        [0, np.sin(i_rad), np.cos(i_rad)]
    ])
    
    R_z_omega = np.array([
        [np.cos(omega_rad), -np.sin(omega_rad), 0],
        [np.sin(omega_rad), np.cos(omega_rad), 0],
        [0, 0, 1]
    ])
    
    rotation_matrix = R_z_Omega @ R_x_i @ R_z_omega
    
    orbit_coords = np.vstack((x_orb, y_orb, z_orb))
    rotated_coords = rotation_matrix @ orbit_coords
    
    return rotated_coords

def main():
    path_to_conf = "venus_orbit_elements.env"  # путь до конфига
    e, T, mu, omega, Omega, i = load_orbit_info(path_to_conf)
    a = (T * T * mu / (4 * np.pi * np.pi)) ** (1/3)  # большая полуось

    orbit_x, orbit_y, orbit_z = orbital_elements_to_cartesian(a, e, i, Omega, omega, T)

    radius_venus = 6051.8  # Радиус Венеры в километрах

    phi, theta = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]  # создание многомерной сетки координат
    x = radius_venus * np.cos(phi) * np.sin(theta)
    y = radius_venus * np.sin(phi) * np.sin(theta)
    z = radius_venus * np.cos(theta)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(x, y, z, cmap=cm.viridis, alpha=0.6)
    ax.set_box_aspect([1,1,1])

    length = (orbit_x ** 2 + orbit_y ** 2 + orbit_z ** 2) ** 0.5
    orbit_x /= length
    orbit_x *= radius_venus

    orbit_y /= length
    orbit_y *= radius_venus

    orbit_z /= length
    orbit_z *= radius_venus

    ax.plot(orbit_x, orbit_y, orbit_z, color='red', label='Траектория The Magellan')

    step = 50
    ax.scatter(orbit_x[::step], orbit_y[::step], orbit_z[::step], color='blue')

    ax.set_xlabel('X (км)')
    ax.set_ylabel('Y (км)')
    ax.set_zlabel('Z (км)')
    ax.set_title('Подспутниковая траектория "The Magellan" вокруг Венеры')
    ax.legend()

    max_val = np.max(radius_venus)
    for axis in 'xyz':
        getattr(ax, f'set_{axis}lim')([-max_val*1.5, max_val*1.5])

    plt.show()

if __name__ == "__main__":
    main()
