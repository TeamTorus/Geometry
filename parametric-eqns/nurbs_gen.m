
%% Nurbs-gen function (call `syms s` first and pass it in as the parameter)
function [x_curve, y_curve, z_curve] = nurbs_gen(s, control_points, weights, to_plot)

    %% Define symbolic control points and weights for 3D NURBS spline
    % control_points = [0 0 0; 5 15 0; 7 5 10; 10 10 5]; % [x, y, z]
    
    % Weights for each control point
    % Keep weights as 1 for a uniform NURBS that passes through the first and last points
    % weights = [1, 1, 1, 1]; 
    
    % Knot vector for NURBS curve (for degree 3 with clamped ends)
    knot_vector = [0 0 0 0 1 1 1 1];
    
    % Number of control points
    n = size(control_points, 1) - 1;
    
    % Degree of the NURBS curve (cubic)
    degree = 3;
    
    % Calculate the symbolic NURBS curve (x, y, z) based on parameter s
    
    x_curve = 0;
    y_curve = 0;
    z_curve = 0;
    
    for i = 1:n+1
        % NURBS basis function for control point i
        N_i = nurbs_basis(i, degree, s, knot_vector);
        
        % Multiply the basis function by the control point and weight
        x_curve = x_curve + N_i * weights(i) * control_points(i, 1);
        y_curve = y_curve + N_i * weights(i) * control_points(i, 2);
        z_curve = z_curve + N_i * weights(i) * control_points(i, 3);
    end
    
    % Divide by the weighted sum of the basis functions
    sum_weights = 0;
    for i = 1:n+1
        sum_weights = sum_weights + nurbs_basis(i, degree, s, knot_vector) * weights(i);
    end
    
    x_curve = x_curve / sum_weights;
    y_curve = y_curve / sum_weights;
    z_curve = z_curve / sum_weights;

    x_curve = simplify(x_curve);
    y_curve = simplify(y_curve);
    z_curve = simplify(z_curve);
    
    %% Plot the NURBS spline
    if to_plot
        s_vals_spline = linspace(0, 1, 100);
    
        x_spline = double(subs(x_curve, s_vals_spline));
        y_spline = double(subs(y_curve, s_vals_spline));
        z_spline = double(subs(z_curve, s_vals_spline));

        figure;
        plot3(x_spline, y_spline, z_spline, 'LineWidth', 2);
        grid on;
        axis equal;
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        title('3D NURBS Spline');
        hold on;
        plot3(control_points(:, 1), control_points(:, 2), control_points(:, 3), '.')
        hold off;
    end
    
    %% Function to compute NURBS basis functions recursively
    function N = nurbs_basis(i, p, u, knot_vector)
        syms N;  % Ensure N is symbolic
        if p == 0
            % Handle symbolic logic using heaviside for basis function
            N = heaviside(u - knot_vector(i)) * heaviside(knot_vector(i+1) - u);
        else
            % First term
            if knot_vector(i+p) == knot_vector(i)
                c1 = 0;
            else
                c1 = ((u - knot_vector(i)) / (knot_vector(i+p) - knot_vector(i))) * nurbs_basis(i, p-1, u, knot_vector);
            end
            
            % Second term
            if knot_vector(i+p+1) == knot_vector(i+1)
                c2 = 0;
            else
                c2 = ((knot_vector(i+p+1) - u) / (knot_vector(i+p+1) - knot_vector(i+1))) * nurbs_basis(i+1, p-1, u, knot_vector);
            end
            
            N = c1 + c2;
        end
    end
end