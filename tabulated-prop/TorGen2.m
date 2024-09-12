airfoilBase = importdata('NACA4415XCT.txt');
xNorm = [airfoilBase(:, 1); flip(airfoilBase(1:length(airfoilBase) - 1, 1))];
camberNorm = [airfoilBase(:, 2); flip(airfoilBase(1:length(airfoilBase) - 1, 2))];
thicknessNorm = [0.5 * airfoilBase(:, 3); -0.5 * flip(airfoilBase(1:length(airfoilBase) - 1, 3))];    

propValues = importdata('torPropValues2.txt');
density = length(propValues(:, 1));
D = 3.6;
R = D/2;
rhoverR = 0.2;
rh = rhoverR * R;
Z = 5; 

roverR = propValues(:, 1);
r = roverR * R;
koverD = propValues(:, 2);
k = koverD * D;
skewDeg = propValues(:, 3);
skewRad = skewDeg * pi /180;
gammaDeg = propValues(:, 4);
boverD = propValues(:, 5);
b = boverD * D;
foverb = propValues(:, 6);
f = foverb .* D;
toverD = propValues(:, 7);
t = toverD .* D;

fileOut = fopen('BladePoints.txt','w');
hold on
axis equal
grid on

guideCurve = zeros(density, 3);

[BCLz, BCLy, BCLx] = pol2cart(skewRad, r, k .* - 1);

for (i = 1:density)
    ax = xNorm * b(i) - b(i)/2;
    ay = camberNorm * f(i) + thicknessNorm * t(i);
    airfoilXY = [ax, ay] * [cosd(-gammaDeg(i)), -sind(-gammaDeg(i)); sind(-gammaDeg(i)), cosd(-gammaDeg(i))]';
    airfoilXYZ = [airfoilXY, zeros(length(airfoilXY), 1)];

    %BCL TNB calculation
    T = [4; 3; 2]; % rng'd
    T = T/norm(T);
    N = [2; 0; -4]; 
    N = N/norm(N);
    B = cross(T, N);
    B = B/norm(B);

    M = [-B, -N, -T];
    Q = airfoilXYZ * M';

    finalAirfoil = Q + [BCLx(i), BCLy(i), BCLz(i)];
    guideCurve(i, :) = finalAirfoil(i);

    plot3(finalAirfoil(:, 1), finalAirfoil(:, 2), finalAirfoil(:, 3), "blue-");
    %plot3(airfoilXYZ(:, 1), airfoilXYZ(:, 2), airfoilXYZ(:, 3), "blue-");
    %plot3(Q(:, 1), Q(:, 2), Q(:, 3), "green-");
    %plot3([0;T(0)], [0;T(1)], [0;T(1)], "red-");
    %plot3([0;N(0)], [0;N(0)], [0;N(0)], "red-");

    fprintf(fileID, "Airfoil%1.0f\n", i);
    fprintf(fileID, "START\t%1f\n", 0);
    for (j = 1:length(finalAirfoil))
        fprintf(fileID, "%9.8f\t%9.8f\t%9.8f\n", finalAirfoil(j, 1), finalAirfoil(j, 2), finalAirfoil(j, 3));
    end
    fprintf(fileID, "END\t%1f\n", 0);
end

plot3(BCLx, BCLy, BCLz, "red-");
xlabel("x or hub");
ylabel("y");
zlabel("z or r");

fprintf(fileID, "GuideCurve%1.0f\n", 1);
fprintf(fileID, "START\t%1f\n", 0);
for (e = 1:density)
    fprintf(fileID, "%9.8f\t%9.8f\t%9.8f\n", guideCurve(e, :));
end
fprintf(fileID, "END\t%1f\n", 0);
fprintf(fileID, "HubRadius %1.f\n", rh);
fprintf(fileID, "Blades %1.f\n", Z);