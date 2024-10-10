import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Define symbolic variable s
s = sp.Symbol('s')

# Nurbs generation function
def nurbs_gen(s, control_points, weights, to_plot=False):
    # Knot vector for NURBS curve (for degree 3 with clamped ends)
    knot_vector = [0, 0, 0, 0, 1, 1, 1, 1]

    # Number of control points
    n = len(control_points) - 1

    # Degree of the NURBS curve (cubic)
    degree = 3

    # Initialize NURBS curve components
    x_curve = 0
    y_curve = 0
    z_curve = 0

    # Calculate NURBS curve based on control points and weights
    for i in range(n + 1):
        # NURBS basis function for control point i
        N_i = nurbs_basis(i, degree, s, knot_vector)
        
        # Multiply basis function by control point and weight
        x_curve += N_i * weights[i] * control_points[i, 0]
        y_curve += N_i * weights[i] * control_points[i, 1]
        z_curve += N_i * weights[i] * control_points[i, 2]

    # Weighted sum of basis functions
    sum_weights = 0
    for i in range(n + 1):
        sum_weights += nurbs_basis(i, degree, s, knot_vector) * weights[i]

    # Divide by weighted sum of basis functions
    x_curve /= sum_weights
    y_curve /= sum_weights
    z_curve /= sum_weights

    # Simplify expressions
    x_curve = sp.simplify(x_curve)
    y_curve = sp.simplify(y_curve)
    z_curve = sp.simplify(z_curve)

    # Plot the NURBS spline if requested
    if to_plot:
        plot_nurbs_curve(x_curve, y_curve, z_curve, control_points)

    return x_curve, y_curve, z_curve


# Function to compute NURBS basis functions recursively
def nurbs_basis(i, p, u, knot_vector):
    if p == 0:
        # Basis function for degree 0 (use heaviside for symbolic logic)
        return sp.Heaviside(u - knot_vector[i]) * sp.Heaviside(knot_vector[i + 1] - u)
    else:
        # First term
        if knot_vector[i + p] == knot_vector[i]:
            c1 = 0
        else:
            c1 = ((u - knot_vector[i]) / (knot_vector[i + p] - knot_vector[i])) * nurbs_basis(i, p - 1, u, knot_vector)

        # Second term
        if knot_vector[i + p + 1] == knot_vector[i + 1]:
            c2 = 0
        else:
            c2 = ((knot_vector[i + p + 1] - u) / (knot_vector[i + p + 1] - knot_vector[i + 1])) * nurbs_basis(i + 1, p - 1, u, knot_vector)

        return c1 + c2


# Function to plot the NURBS spline
def plot_nurbs_curve(x_curve, y_curve, z_curve, control_points):
    s_vals_spline = np.linspace(0, 1, 100)

    # Evaluate the NURBS curve at sampled points
    x_spline = np.array([float(x_curve.subs(s, val)) for val in s_vals_spline])
    y_spline = np.array([float(y_curve.subs(s, val)) for val in s_vals_spline])
    z_spline = np.array([float(z_curve.subs(s, val)) for val in s_vals_spline])

    # 3D plot of the NURBS spline
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x_spline, y_spline, z_spline, label='NURBS Spline', linewidth=2)
    
    # Plot the control points
    ax.plot(control_points[:, 0], control_points[:, 1], control_points[:, 2], 'o', label='Control Points')

    # Set labels and grid
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('3D NURBS Spline')
    plt.legend()
    plt.grid(True)
    plt.show()


# Example usage
if __name__ == "__main__":
    # Define control points as a numpy array (x, y, z)
    control_points = np.array([
        [0, 0, 0],
        [5, 15, 0],
        [7, 5, 10],
        [10, 10, 5]
    ])

    # Define weights for each control point
    weights = [1, 1, 1, 1]  # Uniform weights for this example

    # Generate NURBS curve and plot it
    x_curve, y_curve, z_curve = nurbs_gen(s, control_points, weights, to_plot=True)
