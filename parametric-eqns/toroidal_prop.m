% Symbolic 3D extrusion of a 2D airfoil shape along a 3D curve using the Joukowsky transform

% Clear workspace
clear;
clc;

% Load symbolic library
syms t s

%% Modifiable parameters

num_blades = 3; % Number of blades

% Airfoil Params
thickness = 0.5;
m = 0.02;
p = 0.4;

% Centerline Params
ctrl_point1 = [3 3 10];
ctrl_point2 = [7 6 8];
blade_vector = [8 7];   % offset between the two endpoints

% Angle of Attack
a_AoA = 0;
b_AoA = 0;
c_AoA = 0;
d_AoA = 2 * pi;
e_AoA = 0;

% Scaling Params
a_scX = 0;
b_scX = 0;
c_scX = 0;
d_scX = 1;
e_scX = 0;

a_scY = 0;
b_scY = 0;
c_scY = 0;
d_scY = 1;
e_scY = 0;



%% Define the 2D Shape
% https://en.wikipedia.org/wiki/NACA_airfoil
x_2D = t * heaviside(1 - t) + (2 - t) * heaviside(t - 1);

yt = 5 * thickness * (0.2969 * sqrt(t) - 0.1260 * t - 0.3516 * t^2 + 0.2843 * t^3 - 0.1036 * t^4);
yc = m / p^2 * (2 * p * t - t^2) * heaviside(p - t) + m / (1 - p)^2 * ((1 - 2*p) + 2*p*t - t^2) * heaviside(t - p);
yu = yt + yc;

% domain shift
yl_1 = yc - yt;
yl = subs(yl_1, t, (2-t));

% use activations for upper/lower
y_2D = yu * heaviside(1 - t) + yl * heaviside(t - 1);

%% Define the 3D Curve

% Parameters for the 3D curve (a helix in this case)
% R = 10; % Radius of the helix
% pitch = 2; % Vertical distance per revolution
% n_turns = 3; % Number of turns

% Parametric equations for the 3D curve 
% using these parabolic eqns use domain of like -2 to 2
% x_curve = R * cos(s);
% y_curve = R * sin(s);
% z_curve = -1 * s ^ 2;
control_points = [0 0 0; ctrl_point1; ctrl_point2; blade_vector(1) blade_vector(2) 0]; % [x, y, z]
[x_curve, y_curve, z_curve] = nurbs_gen(s, control_points, [1,1,1,1], false)

%% Define the Frenet-Serret Frame

% Derivatives to get tangent vector (T)
dx_ds = diff(x_curve, s);
dy_ds = diff(y_curve, s);
dz_ds = diff(z_curve, s);
T = [dx_ds, dy_ds, dz_ds];
T = T / norm(T);

% Derivative of tangent to get normal vector (N)
dT_ds = diff(T, s);
N = dT_ds / norm(dT_ds);

% Binormal vector (B) is the cross product of T and N
B = cross(T, N);

%% Define Rotation Angle based on curve parameter s
AoA = a_AoA * s^4 + b_AoA * s^3 + c_AoA * s^2 + d_AoA * s + e_AoA;

% Define the rotation matrix (counterclockwise in the local x-y plane)
rotation_matrix = [cos(AoA), -sin(AoA); sin(AoA), cos(AoA)];

% Rotate the 2D shape symbolically
XY_rotated = rotation_matrix * [x_2D; y_2D];  % Apply rotation
X_rotated = XY_rotated(1);
Y_rotated = XY_rotated(2);

%% Define symbolic scaling/stretching functions (parametric by s)
scale_x = a_scX * s^4 + b_scX * s^3 + c_scX * s^2 + d_scX * s + e_scX;  % Parametric scaling for x
scale_y = a_scY * s^4 + b_scY * s^3 + c_scY * s^2 + d_scY * s + e_scY;  % Parametric scaling for y

% Apply scaling in the local frame (after rotation)
X_rotated_scaled = X_rotated * scale_x;
Y_rotated_scaled = Y_rotated * scale_y;

%% Transform the 2D shape Along the Curve

% Define the parametric position along the curve
C = [x_curve, y_curve, z_curve];

% Express the final 3D coordinates of the airfoil after extrusion
X_final = C(1) + X_rotated_scaled * N(1) + Y_rotated_scaled * B(1);
Y_final = C(2) + X_rotated_scaled * N(2) + Y_rotated_scaled * B(2);
Z_final = C(3) + X_rotated_scaled * N(3) + Y_rotated_scaled * B(3);

% Convert symbolic expressions to functions for numerical evaluation
X_func = matlabFunction(X_final, 'Vars', [s, t]);
Y_func = matlabFunction(Y_final, 'Vars', [s, t]);
Z_func = matlabFunction(Z_final, 'Vars', [s, t]);

%% Numerical Evaluation and Plotting

% Discretize the parameters
s_vals = linspace(0, 1);  % Values for curve parameter (s)
t_vals = linspace(0, 2, 1000);  % Values for shape parameter (t)

% Evaluate the symbolic expressions numerically
[X_mesh, T_mesh] = meshgrid(s_vals, t_vals);
X_vals = X_func(X_mesh, T_mesh);
Y_vals = Y_func(X_mesh, T_mesh);
Z_vals = Z_func(X_mesh, T_mesh);

% Plot the extruded prop in 3D
figure;
surf(X_vals, Y_vals, Z_vals, 'EdgeColor', 'none');
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Extruded Toroidal Prop Blade');
camlight; lighting phong;
