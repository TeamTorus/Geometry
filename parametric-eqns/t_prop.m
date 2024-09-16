syms x y z real

%% Define two cylinders

% Cylinder 1 (aligned along the z-axis)
r1 = 1; % radius of cylinder 1
D1 = sqrt(x^2 + y^2) - r1; % Signed distance function for cylinder 1

% Cylinder 2 (aligned along the y-axis)
r2 = 1; % radius of cylinder 2
D2 = sqrt(x^2 + z^2) - r2; % Signed distance function for cylinder 2

% Union of the two cylinders using the minimum distance
D_union = min(D1, D2)

%% Visualize the result
% Generate a grid of points to evaluate the union distance function
[x_vals, y_vals, z_vals] = meshgrid(linspace(-2, 2, 50), linspace(-2, 2, 50), linspace(-2, 2, 50));

% Evaluate the union distance function at each point
D_union_vals = double(subs(D_union, {x, y, z}, {x_vals, y_vals, z_vals}));

% Create an isosurface plot for the union of the two cylinders
figure;
isosurface(x_vals, y_vals, z_vals, D_union_vals, 0); % Plot the zero level set
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Union of Two Intersecting Cylinders');
camlight; lighting phong;
