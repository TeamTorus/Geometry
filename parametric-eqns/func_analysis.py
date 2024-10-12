import inspect
import numpy as np
import random
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# define function types to randomly generate
def random_sinusoidal(x):
    amplitude = random.uniform(0.5, 2)
    frequency = random.uniform(1, 3)
    phase = random.uniform(0, np.pi)
    return amplitude * np.sin(2 * np.pi * frequency * x + phase)

def random_polynomial(x):
    degree = random.randint(1, 8)  # Random degree polynomial between 2 and 5
    coefficients = [random.uniform(-2, 2) for _ in range(degree + 1)]
    return sum(c * x**i for i, c in enumerate(coefficients))

# can approximate a square function, or a sharper transition @ blade edge
def random_logistic(x):
    L = random.uniform(1, 50)  # Carrying capacity
    k = random.uniform(8, 100)  # Growth rate
    x0 = random.uniform(0.1, 0.9)  # Midpoint
    return L / (1 + np.exp(-k * (x - x0)))

def random_distribution(x):

    # downsample function to 2-8 interp points to help it fit (otherwise too chaotic)
    n_points = random.randint(3, 8)
    x_points = np.linspace(0, 1, n_points)
    y_points = np.random.uniform(-1, 1, n_points)
    return np.interp(x, x_points, y_points)

func_generators = {
    "sinusoidal": random_sinusoidal,
    "polynomial": random_polynomial,
    "logistic": random_logistic,
    "distribution": random_distribution
}

# Generate random function and data
def generate_random_function(x):
    func_type = random.choice(list(func_generators.keys()))
    func = func_generators[func_type]
    return func(x)

# Polynomial function to fit
def poly4(x, a, b, c, d, e):
    return a * x**4 + b * x**3 + c * x**2 + d * x + e

def poly5(x, a, b, c, d, e, f):
    return a * x**5 + b * x**4 + c * x**3 + d * x**2 + e * x + f

def poly6(x, a, b, c, d, e, f, g):
    return a * x**6 + b * x**5 + c * x**4 + d * x**3 + e * x**2 + f * x + g

# N=2 Fourier Series function to fit
def fourier2(x, A1, B1, A2, B2, A0):
    return A0 + A1 * np.cos(2 * np.pi * x) + B1 * np.sin(2 * np.pi * x) + A2 * np.cos(4 * np.pi * x) + B2 * np.sin(4 * np.pi * x)

def fourier3(x, A1, B1, A2, B2, A3, B3, A0):
    return A0 + A1 * np.cos(2 * np.pi * x) + B1 * np.sin(2 * np.pi * x) + A2 * np.cos(4 * np.pi * x) + B2 * np.sin(4 * np.pi * x) + A3 * np.cos(6 * np.pi * x) + B3 * np.sin(6 * np.pi * x)

def sinusoidal(x, a, b, c, d):
    return a * np.sin(b * x + c) + d

fitting_functions = {
    "polynomial4": poly4,
    "polynomial6": poly6,
    "fourier2": fourier2,
    "fourier3": fourier3,
    "sinusoidal": sinusoidal
}

# fit polynomial and Fourier series, using MSE metric (could be optimized or regularized but idrk)
def fit_models(x_data, y_data):

    fit_outputs = {}

    # Fit the funcs
    for func_name, func in fitting_functions.items():
        try:
            params, _ = curve_fit(func, x_data, y_data, maxfev=10000)
        except RuntimeError:        
            # If the fit fails, set all params to 0
            params = np.zeros(len(inspect.signature(func).parameters) - 1)
            params[len(params) - 1] = np.mean(y_data)       # mean at data to avoid crazy overfitted errors
        
        # Get residuals
        residuals = y_data - func(x_data, *params)

        # Normalize residuals against max possible residual
        residuals /= (np.max(y_data) - np.min(y_data))
        error = np.sum(residuals**2)

        fit_outputs[func_name] = (error, params)

    return fit_outputs

n_samples = 10000
n_points = 100  # Number of test data points

# dictionary comprehension to store errors
total_errors = {func_name: [] for func_name in fitting_functions.keys()}

# x is [0, 1]
x_data = np.linspace(0, 1, n_points)

# iterate
for _ in range(n_samples):
    # Generate random function data
    y_data = generate_random_function(x_data)
    
    # Fit both models and record errors
    fit = fit_models(x_data, y_data)
    
    for func_name, (error, _) in fit.items():
        total_errors[func_name].append(error)

    # plot the test data and the fitted models
    if n_samples == 1:

        plt.figure()
        plt.plot(x_data, y_data, label="Test Data")

        for func_name, (error, params) in fit.items():
            print(func_name, "Parameters: ", params)

            plt.plot(x_data, fitting_functions[func_name](x_data, *params), label=func_name)

        plt.legend()
        plt.show()

for func_name, errors in total_errors.items():
    print(func_name, "Mean Error:", np.mean(errors))

