% Symbolic 3D extrusion of a 2D airfoil shape along a 3D curve using the Joukowsky transform

% Clear workspace
clear;
clc;

% Load symbolic library
syms t s

s_domain = [0, 1]; % Domain for the curve parameter s
t_domain = [0, 2]; % Domain for the shape parameter t
s_resolution = 100; % Resolution for discretizing the curve parameter s
t_resolution = 1000; % Resolution for discretizing the shape parameter t

%% Modifiable parameters

hub_radius = 5; % Radius of the cylindrical hub
hub_length = 20; % Length of the cylindrical hub
num_blades = 3; % Number of blades

% Airfoil Params
thickness = 0.5;
m = 0.02;
p = 0.4;

% Centerline Params
loc_ctrl_point2 = [3 3 10];
loc_ctrl_point3 = [7 6 8];
blade_vector = [8 8];   % offset between the two endpoints

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
e_scX = 2;

a_scY = 0;
b_scY = 0;
c_scY = 0;
d_scY = 1;
e_scY = 2;



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

% pick first point
ctrl_point1 = [hub_radius, 0, hub_length/2 - 1];

% for last control point, offset by blade_vector
% blade_vector(1) represents a movement along the circumference of the hub
% blade_vector(2) represents a movement along the z-axis
disp_theta = blade_vector(1) / (2 * pi * hub_radius) * 2 * pi;
ctrl_point4 = [hub_radius * cos(disp_theta), hub_radius * sin(disp_theta), -1 * blade_vector(2)];

% Define the control points for the curve (converting to cylindrical and translating)
disp_theta2 = loc_ctrl_point2(1) / (2 * pi * hub_radius) * 2 * pi;
loc2_radius = hub_radius + loc_ctrl_point2(3);
ctrl_point2 = [loc2_radius * cos(disp_theta2), loc2_radius * sin(disp_theta2), -1 * loc_ctrl_point2(2)];

disp_theta3 = loc_ctrl_point3(1) / (2 * pi * hub_radius) * 2 * pi;
loc3_radius = hub_radius + loc_ctrl_point3(3);
ctrl_point3 = [loc3_radius * cos(disp_theta3), loc3_radius * sin(disp_theta3), -1 * loc_ctrl_point3(2)];

control_points = [ctrl_point1; ctrl_point2; ctrl_point3; ctrl_point4]; % [x, y, z]
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

%% Assemble the 3D Propeller Blade

% convert to cylindrical coordinates
R = sqrt(X_final^2 + Y_final^2);
Theta = atan2(Y_final, X_final);

% create num_blades - 1 additional blades, rotating around the hub

X_prop = zeros(t_resolution, s_resolution * num_blades);
Y_prop = zeros(t_resolution, s_resolution * num_blades);
Z_prop = zeros(t_resolution, s_resolution * num_blades);


for i = 0:num_blades - 1

    % rotate the blade by 2pi/num_blades
    Loc_Theta = Theta + i * 2 * pi / num_blades;
    
    % convert back
    X_rotated = R * cos(Loc_Theta);
    Y_rotated = R * sin(Loc_Theta);

    %% Numerical Evaluation and Plotting

    % Convert symbolic expressions to functions for numerical evaluation
    X_func = matlabFunction(X_rotated, 'Vars', [s, t]);
    Y_func = matlabFunction(Y_rotated, 'Vars', [s, t]);
    Z_func = matlabFunction(Z_final, 'Vars', [s, t]);

    % Discretize the parameters for plotting
    s_vals = linspace(s_domain(1), s_domain(2), s_resolution);  % Values for curve parameter (s)
    t_vals = linspace(t_domain(1), t_domain(2), t_resolution);  % Values for shape parameter (t)

    % Evaluate the symbolic expressions numerically
    [X_mesh, T_mesh] = meshgrid(s_vals, t_vals);
    X_vals = X_func(X_mesh, T_mesh);
    Y_vals = Y_func(X_mesh, T_mesh);
    Z_vals = Z_func(X_mesh, T_mesh);

    X_prop(:, i * s_resolution + 1:(i + 1) * s_resolution) = X_vals;
    Y_prop(:, i * s_resolution + 1:(i + 1) * s_resolution) = Y_vals;
    Z_prop(:, i * s_resolution + 1:(i + 1) * s_resolution) = Z_vals;

end

% Generate Cylinder Hub (parametric)
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

%% Perform Union: Remove points inside the hub

% Create a mask for points inside the hub
distance_from_origin = sqrt(X_prop.^2 + Y_prop.^2);
inside_hub_mask = (distance_from_origin < hub_radius);

% Set the points inside the hub to NaN (removes them)
X_prop(inside_hub_mask) = NaN;
Y_prop(inside_hub_mask) = NaN;
Z_prop(inside_hub_mask) = NaN;

% Plot the extruded prop in 3D
surf(X_prop, Y_prop, Z_prop, 'EdgeColor', 'none');
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Combined Propeller Blades and Hub');
camlight; lighting phong;