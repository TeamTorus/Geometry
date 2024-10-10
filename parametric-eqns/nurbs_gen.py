import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Define symbolic variable s
s = sp.Symbol('s')

# NURBS generation function (translated from MATLAB)
def nurbs_gen(s, control_points, weights, to_plot=False):
    # cubic nurbs
    degree = 3

    # Knot vector for NURBS curve (for degree 3 with clamped ends)
    knot_vector = np.array([0, 0, 0, 0, 1, 1, 1, 1])

    n = len(control_points) - 1

    def nurbs_basis(i, p, u, knot_vector):
        if p == 0:
            return sp.Piecewise((1, (u >= knot_vector[i]) & (u < knot_vector[i+1])), (0, True))
        else:
            denom1 = knot_vector[i+p] - knot_vector[i]
            denom2 = knot_vector[i+p+1] - knot_vector[i+1]
            term1 = 0
            term2 = 0
            if denom1 != 0:
                term1 = ((u - knot_vector[i]) / denom1) * nurbs_basis(i, p-1, u, knot_vector)
            if denom2 != 0:
                term2 = ((knot_vector[i+p+1] - u) / denom2) * nurbs_basis(i+1, p-1, u, knot_vector)
            return term1 + term2

    x_curve = 0
    y_curve = 0
    z_curve = 0
    sum_weights = 0

    for i in range(n + 1):
        N_i = nurbs_basis(i, degree, s, knot_vector)
        wN_i = N_i * weights[i]
        x_curve += wN_i * control_points[i, 0]
        y_curve += wN_i * control_points[i, 1]
        z_curve += wN_i * control_points[i, 2]
        sum_weights += wN_i

    x_curve = x_curve / sum_weights
    y_curve = y_curve / sum_weights
    z_curve = z_curve / sum_weights

    x_curve = sp.simplify(x_curve)
    y_curve = sp.simplify(y_curve)
    z_curve = sp.simplify(z_curve)

    if to_plot:
        plot_nurbs_curve(x_curve, y_curve, z_curve, control_points)

    return x_curve, y_curve, z_curve


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


# usage
if __name__ == "__main__":
    # control points as a numpy array (x, y, z)
    control_points = np.array([
        [0, 0, 0],
        [5, 15, 0],
        [7, 5, 10],
        [10, 10, 5]
    ])

    weights = [1, 1, 1, 1]  # uniform weights

    # generate NURBS curve and plot it
    x_curve, y_curve, z_curve = nurbs_gen(s, control_points, weights, to_plot=True)
