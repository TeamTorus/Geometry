propValues = importdata('0000BladeParameters.txt');
fileID = fopen('TorBladePoints.txt','w');

D = 30.48;
LoverD = 0.3;
L = LoverD * D;
R = D/2;
rhoverR = 0.2;
rh = rhoverR * R;  
Z = 3;  
loverL = propValues(:, 1);
l = loverL * L;
rloverR = propValues(:, 2);
rl = rloverR * R;
boverD = propValues(:, 3);
b = boverD * D;
PoverD = propValues(:, 4);
P = PoverD * D;
beta = atan(P ./ (2 * pi * rl));
toverb = propValues(:, 5);
t = toverb .* b;
foverb = propValues(:, 6);
f = foverb .* b;
thetaSDeg = propValues(:, 7);
thetaS = thetaSDeg * pi / 180;
xloverD = propValues(:, 8);
xl = xloverD * D;
phiDeg = propValues(:, 9);
phi = phiDeg  * pi / 180;
psiDeg = propValues(:, 10);
psi = psiDeg  * pi / 180;
alphaDeg = propValues(:, 11);
alpha = alphaDeg * pi / 180;

density = length(rl);
airfoilPoints = 53;
guideCurves = zeros(density, 3, airfoilPoints);

% hold on
% axis equal
for (i = 1:density)
    [s, y] = NACA66TMBModBuilder(b(i), t(i), f(i));

    r = rl(i) - (-b(i)/2 + s) * sin(alpha(i)) - y * sin(psi(i));
    c = b(i)/2 - ((r * thetaS(i)) / (cos(beta(i)) * cos(alpha(i))));
    theta = phi(i) + (1./r) .* (((-c + s) * cos(alpha(i))) * cos(beta(i)) + (y * cos(psi(i)) * sin(beta(i))));
    x = l(i) + xl(i) + ((-c + s) * cos(alpha(i))) * sin(beta(i)) - (y * cos(psi(i))) * cos(beta(i));

    xPlot = x;
    yPlot = r .* cos(theta);
    zPlot = r .* sin(theta); 
    % if (i == 8 || i == 16)
    %     plot3(xPlot, yPlot, zPlot, 'blue-');
    % elseif (i == 12)
    %     plot3(xPlot, yPlot, zPlot, 'green-');
    % else
    %     plot3(xPlot, yPlot, zPlot, 'red-');
    % end

    guideCurves(i, 1, :) = xPlot();
    guideCurves(i, 2, :) = yPlot();
    guideCurves(i, 3, :) = zPlot();
    fprintf(fileID, "Airfoil%1.0f\n", i);
    fprintf(fileID, "START\t%1f\n", 0);
    for (j = 1:length(s))
        fprintf(fileID, "%9.8f\t%9.8f\t%9.8f\n", xPlot(j), yPlot(j), zPlot(j));
    end
    fprintf(fileID, "END\t%1f\n", 0);
end

for (g = 1:airfoilPoints)
    fprintf(fileID, "GuideCurve%1.0f\n", g);
    fprintf(fileID, "START\t%1f\n", 0);
    for (e = 1:density)
        fprintf(fileID, "%9.8f\t%9.8f\t%9.8f\n", guideCurves(e, :, g));
    end
    fprintf(fileID, "END\t%1f\n", 0);
end

fprintf(fileID, "HubRadius %2.5f\n", rh);
fprintf(fileID, "Blades %1.f\n", Z);