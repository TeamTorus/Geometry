% Clear workspace
clear;
clc;

% Load symbolic library
syms t s

%% Modifiable parameters

num_blades = 3; % Number of blades
hub_radius = 10; % Radius of the cylindrical hub
hub_length = 20; % Length of the cylindrical hub

% centerline params
blade_vector = [8 7]; % Offset between the two endpoints of blade (in section coordinates)
loc_ctrl_point2 = [3 3 10];
loc_ctrl_point3 = [7 6 8];

% Airfoil Params
thickness = 0.5;
m = 0.02;
p = 0.4;

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

%% Generate Cylinder Hub
theta_hub = linspace(0, 2*pi, 100);
z_hub = linspace(-hub_length/2, hub_length/2, 50);
[TH, ZH] = meshgrid(theta_hub, z_hub);

X_hub = hub_radius * cos(TH);
Y_hub = hub_radius * sin(TH);
Z_hub = ZH;

% Plot hub
figure;
surf(X_hub, Y_hub, Z_hub, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
axis equal;

%% Generate Blades Around the Hub

for k = 0:num_blades-1
    % Rotation angle for current blade section (evenly spaced around the hub)
    rotation_angle = 2 * pi * k / num_blades;
    
    % Top left corner of the current section (start of spline)
    ctrl_point1 = hub_radius * [cos(rotation_angle), sin(rotation_angle), hub_length/2];
    
    % Offset the last control point by blade_vector relative to the section
    ctrl_point4 = ctrl_point1 + [blade_vector(1)*cos(rotation_angle), blade_vector(1)*sin(rotation_angle), -blade_vector(2)];
    
    control_points = [ctrl_point1; (ctrl_point1 + [3 3 3]); (ctrl_point1 + [5 5 5]); ctrl_point4]; % [x, y, z]

    % NURBS spline for the current blade
    [x_curve, y_curve, z_curve] = nurbs_gen(s, control_points, [1,1,1,1], false);
    
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
    
    % Define Rotation Angle based on curve parameter s (angle of attack)
    AoA = a_AoA * s^4 + b_AoA * s^3 + c_AoA * s^2 + d_AoA * s + e_AoA;
    
    % Rotation matrix for blade orientation
    rotation_matrix = [cos(AoA), -sin(AoA); sin(AoA), cos(AoA)];
    
    % Define 2D Airfoil shape
    x_2D = t * heaviside(1 - t) + (2 - t) * heaviside(t - 1);
    yt = 5 * thickness * (0.2969 * sqrt(t) - 0.1260 * t - 0.3516 * t^2 + 0.2843 * t^3 - 0.1036 * t^4);
    yc = m / p^2 * (2 * p * t - t^2) * heaviside(p - t) + m / (1 - p)^2 * ((1 - 2*p) + 2*p*t - t^2) * heaviside(t - p);
    yu = yt + yc;
    yl_1 = yc - yt;
    yl = subs(yl_1, t, (2-t));
    y_2D = yu * heaviside(1 - t) + yl * heaviside(t - 1);
    
    % Rotate the 2D shape symbolically
    XY_rotated = rotation_matrix * [x_2D; y_2D];  % Apply rotation
    X_rotated = XY_rotated(1);
    Y_rotated = XY_rotated(2);
    
    % Define symbolic scaling/stretching functions (parametric by s)
    scale_x = a_scX * s^4 + b_scX * s^3 + c_scX * s^2 + d_scX * s + e_scX;
    scale_y = a_scY * s^4 + b_scY * s^3 + c_scY * s^2 + d_scY * s + e_scY;
    
    % Apply scaling in the local frame (after rotation)
    X_rotated_scaled = X_rotated * scale_x;
    Y_rotated_scaled = Y_rotated * scale_y;
    
    % Express the final 3D coordinates of the airfoil after extrusion
    X_final = x_curve + X_rotated_scaled * N(1) + Y_rotated_scaled * B(1);
    Y_final = y_curve + X_rotated_scaled * N(2) + Y_rotated_scaled * B(2);
    Z_final = z_curve + X_rotated_scaled * N(3) + Y_rotated_scaled * B(3);
    
    % Convert symbolic expressions to functions for numerical evaluation
    X_func = matlabFunction(X_final, 'Vars', [s, t]);
    Y_func = matlabFunction(Y_final, 'Vars', [s, t]);
    Z_func = matlabFunction(Z_final, 'Vars', [s, t]);
    
    % Discretize parameters for numerical evaluation
    s_vals = linspace(0, 1);
    t_vals = linspace(0, 2, 1000);
    
    [X_mesh, T_mesh] = meshgrid(s_vals, t_vals);
    X_vals = X_func(X_mesh, T_mesh);
    Y_vals = Y_func(X_mesh, T_mesh);
    Z_vals = Z_func(X_mesh, T_mesh);
    
    % Plot the current blade
    surf(X_vals, Y_vals, Z_vals, 'EdgeColor', 'none');
end

hold off;
title('Full Toroidal Propeller with Multiple Blades');
xlabel('X'); ylabel('Y'); zlabel('Z');
camlight; lighting phong;
