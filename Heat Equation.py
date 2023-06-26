import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def solve_heat_equation(length, max_time, num_steps_x, num_steps_t, thermal_diffusivity):
    """
    Solve the 1D heat equation using the explicit method and visualize the solution.
    
    Parameters:
        length (float): Length of the rod.
        max_time (float): Maximum time.
        num_steps_x (int): Number of steps in the x-direction.
        num_steps_t (int): Number of steps in the t-direction.
        thermal_diffusivity (float): Thermal diffusivity.
    """
    # Create the grid
    x = np.linspace(0, length, num_steps_x + 1)
    t = np.linspace(0, max_time, num_steps_t + 1)
    dx = x[1] - x[0]
    dt = t[1] - t[0]
    u = np.zeros((num_steps_t + 1, num_steps_x + 1))

    # Set the initial condition
    u[0, :] = np.sin(np.pi * x / length)

    # Set the boundary conditions
    u[:, 0] = 0
    u[:, -1] = 0

    # Solve the heat equation using the explicit method
    alpha = thermal_diffusivity
    for j in range(num_steps_t):
        for i in range(1, num_steps_x):
            u[j + 1, i] = u[j, i] + alpha * dt / dx**2 * (u[j, i - 1] - 2 * u[j, i] + u[j, i + 1])

    # Plot the solution
    X, T = np.meshgrid(x, t)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, T, u)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('u')
    plt.show()

# Set the parameters
length = 1
max_time = 0.1
num_steps_x = 100
num_steps_t = 100
thermal_diffusivity = 0.01

# Solve and visualize the heat equation
solve_heat_equation(length, max_time, num_steps_x, num_steps_t, thermal_diffusivity)