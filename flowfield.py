from cmath import sqrt

import numpy as np
import matplotlib.pyplot as plt


# Given equations for the flow rates and stream function contributions
def Q0(a):
    return ((1 - 4 * a ** 2) ** (5 / 2)) / (1 + 2 * a ** 2)


def Q2(a):
    return -12 * (1 - 4 * a ** 2) ** (7 / 2) * np.pi * a ** 2 / (5 * (1 + 2 * a ** 2) ** 2)


def Q4(a):
    return (144 * (1 - 4 * a ** 2) ** (9 / 2) * np.pi * a ** 4 / (25 * (1 + 2 * a ** 2) ** 3) -
            (1 - 4 * a ** 2) ** (5 / 2) * np.pi ** 4 / (175 * (1 + 2 * a ** 2) ** 2) *
            (7648 * a ** 6 - 7680 * a ** 4 + 2406 * a ** 2) +
            214 * (1 - 4 * a ** 2) ** (5 / 2) / 241)


def Psi_0(xi):
    return Q0(a) / 4 * (3 * xi - xi ** 3)


def Psi_2(eta, xi, a):
    # Using dh_deta from the previous function
    h_eta = 1 / 2 + a * np.sin(2 * np.pi * eta)
    dh_deta = a * 2 * np.pi * np.cos(2 * np.pi * eta)
    return Q2(a) / Q0(a) * Psi_0(xi) + 3 * Q0(a) / 40 * (4 * (dh_deta ** 2) - (h_eta * dh_deta)) * xi * ((xi ** 2 - 1) ** 2)


# Update stream_function to include the new terms
def stream_function(eta, xi, epsilon, a):
    # Channel height function and its derivatives
    eta = eta * epsilon
    h_eta = 1 / 2 + a * np.sin(2 * np.pi * eta)
    xi = xi / h_eta
    dh_deta = a * 2 * np.pi * np.cos(2 * np.pi * eta)
    d2h_deta2 = -4 * np.pi ** 2 * a * np.sin(2 * np.pi * eta)
    d3h_deta3 = -8 * np.pi ** 3 * a * np.cos(2 * np.pi * eta)
    d4h_deta4 = 16 * np.pi ** 4 * a * np.sin(2 * np.pi * eta)

    # Zeroth, second, and fourth-order stream function
    psi_0 = Psi_0(xi)
    psi_2 = Psi_2(eta, xi, a)
    psi_4 = ((Q4(a) / Q0(a)) - (Q2(a) / Q0(a)) ** 2) * psi_0 + (Q2(a) / Q0(a)) * psi_2 - (Q0(a) / 5600 * ((408 + 1800 * (xi ** 2)) * (dh_deta ** 4) +
                             (684 - 1800 * (xi ** 2)) * h_eta * (dh_deta ** 2) * d2h_deta2 -
                             (270 - 180 * (xi ** 2)) * (h_eta ** 2) * (d2h_deta2 ** 2) -
                             (248 - 240 * (xi ** 2)) * dh_deta * d3h_deta3 * (h_eta ** 2) +
                             (19 - 15 * (xi ** 2)) * (h_eta ** 3) * d4h_deta4)) * xi * (((xi ** 2) - 1) ** 2)

    return psi_0 + psi_2 + psi_4


# Define the aspect ratio a and fluctuation ratio epsilon
a = 0.4
epsilon = 0.4

# Create a meshgrid for the eta and xi coordinates
x = np.linspace(-2, 5, 30)
y = np.linspace(-1, 1, 30)
eta, xi = np.meshgrid(x, y)

# Calculate the stream function values
Psi = stream_function(eta, xi, a, epsilon)

# Calculate the velocity field from the stream function
Vx = np.gradient(Psi, axis=1)
Vy = -np.gradient(Psi, axis=0)

# Plot the streamlines
plt.figure(figsize=(8, 8))
plt.streamplot(x, y, Vy, Vx, color='black', density=2)
plt.title('Streamlines for given aspect ratio (a) and fluctuation ratio (epsilon)')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

"""

# Aspect ratio and fluctuation ratio
A = 0.2  # Aspect ratio (height to length scale of waves)
F = 0.2  # Fluctuation ratio (amplitude to mean height of boundary)

# Define the grid in x and y dimensions
x = np.linspace(-2, 5, 20)
y = np.linspace(-1, 1, 20)
X, Y = np.meshgrid(x, y)

# Parameters for the wavy wall
L = 7.0  # Length of the domain in the x-direction
h_bar = A * L / 2.0
h_apostrophe = F * h_bar * 2.0
h_x = h_bar + h_apostrophe * np.sin(2 * np.pi * X / L)

# Stream function


# Calculate velocity components from the stream function
u_new = np.gradient(h_x, axis=1) / np.gradient(y, axis=0)  # derivative with respect to y
v_new = np.gradient(h_x, axis=1) / np.gradient(x, axis=0)  # derivative with respect to x

# Plotting the vector field
plt.figure(figsize=(10, 10))
plt.quiver(X, Y, u_new, v_new, scale=20)
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Vector Field of Solute Dispersion with Ratios')
plt.grid()
plt.show()


def compute_velocity_with_ratios(x, y, aspect_ratio, fluctuation_ratio, flow_speed=1.0):
    Compute the velocity of solute particles at a given point (x, y) with aspect ratio and fluctuation ratio.

    :param x: x-coordinate
    :param y: y-coordinate
    :param aspect_ratio: Aspect ratio (2 * avg channel height / L)
    :param fluctuation_ratio: Fluctuation ratio (amplitude of fluctuations / 2 * avg channel height)
    :param flow_speed: Speed of advection (default 1.0)
    :return: Velocity components (vx, vy)

    # Assuming some values for demonstration

    avg_channel_height = aspect_ratio * 5.5 / 2.0
    L = 5.0  # Length of the domain in the x-direction
    if aspect_ratio == 0: amplitude_of_fluctuations = 0
    else: amplitude_of_fluctuations = fluctuation_ratio * 2 * avg_channel_height
    height = avg_channel_height + amplitude_of_fluctuations * np.sin(2 * np.pi * x / L)
    Psi = flow_speed * y * (height - y) / (height ** 2)
    vx = np.ones_like(X) * flow_speed
    vy = -np.gradient(Psi, axis=1) / np.gradient(x, axis=0)
    # Advection: Uniform flow along the x-axis, modified by aspect ratio
    # Diffusion: Velocity fluctuations based on fluctuation ratio
    # vy = fluctuation_ratio * np.sin(np.pi * y / aspect_ratio)

    return vx, vy


# Calculating aspect ratio and fluctuation ratio
# aspect_ratio = 2 * avg_channel_height / L
# fluctuation_ratio = amplitude_of_fluctuations / (2 * avg_channel_height)

aspect_ratio = 0.2
fluctuation_ratio = 0.2

# Define a grid
x = np.linspace(-1.5, 4, 20)
y = np.linspace(-1, 1, 20)
X, Y = np.meshgrid(x, y)

# Compute the velocity field
VX, VY = np.vectorize(compute_velocity_with_ratios)(X, Y, aspect_ratio, fluctuation_ratio)

# Plotting the vector field
plt.figure(figsize=(10, 10))
plt.quiver(X, Y, VX, VY, scale=30)
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Vector Field of Solute Dispersion with Ratios')
plt.grid()
plt.show()
"""
