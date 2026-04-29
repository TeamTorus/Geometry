D = 30.48;
R = D/2;
rhoverR = 0.2;
rh = rhoverR * R;
roverR = [.2, .3, .4, .5, .6, .7, .8, .9, .95, 1];
r = roverR * R;
Z = 3;
density = length(r);

% First 4119 Paper: https://apps.dtic.mil/sti/pdfs/AD0841845.pdf
%coverD = [0.32, 0.3635, 0.4048, 0.4392, 0.461, 0.4622, 0.4347, 0.3613, 0.2775, .01];
%c = coverD * D;
%toverD = [0.0878, 0.0564, 0.0478, 0.0396, 0.032, 0.025, 0.0182, 0.012, 0.009, .0003];
%t = toverD * D;
%foverC = [0.0224, 0.0232, 0.0233, 0.0218, 0.0207, 0.0200, 0.0197, 0.0182, 0.0167, 0.0167];
%f = foverC .* c;
%PoverD = [1.105, 1.102, 1.098, 1.093, 1.088, 1.084, 1.081, 1.078, 1.075, 1.075];
%P = PoverD * D;
%alphaRad = atan(P ./ (2 * pi * r));
%alphaDeg = alphaRad * 180 / pi;
%alphat = [0.952, 0.912, 0.79, 0.612, 0.45, 0.324, 0.23, 0.18, 0.15, 0.15];

% Second 4119 Paper: https://apps.dtic.mil/sti/tr/pdf/ADA222492.pdf
coverD = [0.32, 0.3635, 0.4048, 0.4392, 0.461, 0.4622, 0.4347, 0.3613, 0.2775, .01];
c = coverD * D;
toverC = [0.2055, 0.1552, 0.118, 0.09016, 0.0696, 0.05418, 0.04206, 0.03321, 0.03228, 0.0316];
t = toverC .* c;
foverC = [0.01429, 0.02318, 0.02303, 0.02182, 0.02072, 0.0200, 0.01967, 0.01817, 0.01631, 0.01175];
f = foverC .* c;
PoverD = [1.105, 1.102, 1.098, 1.093, 1.088, 1.084, 1.081, 1.079, 1.077, 1.075];
P = PoverD * D;
alphaRad = atan(P ./ (2 * pi * r));
alphaDeg = alphaRad * 180 / pi;

fileID = fopen('4119BladePoints.txt','w');
guideCurveCount = 53;
guideCurves = zeros(density, 3, guideCurveCount);

for (i = 1:density)
    [airfoilx, airfoily] = NACA66TMBModBuilder(c(i), t(i), f(i));
    airfoilx = airfoilx - c(i)/2;
    airfoilx = airfoilx .* -1;
    airfoil = [airfoilx, airfoily];
    [x, y, z] = Wrapping3D(airfoil, r(i), 1, -1 * alphaDeg(i), 0, 0);
    for (k = 1:guideCurveCount)
        guideCurves(i, :, k) = [x(k), y(k), z(k)];
    end

    fprintf(fileID, "Airfoil%1.0f\n", i);
    fprintf(fileID, "START\t%1f\n", 0);
    for (j = 1:length(x))
        fprintf(fileID, "%9.8f\t%9.8f\t%9.8f\n", x(j), y(j), z(j));
    end
    fprintf(fileID, "END\t%1f\n", 0); 
end

fprintf(fileID, "TangencyLine%1.0f\n", 1);
fprintf(fileID, "START\t%1f\n", 0); 
fprintf(fileID, "%9.8f\t%9.8f\t%9.8f\n", guideCurves(density, :, 1));
fprintf(fileID, "%9.8f\t%9.8f\t%9.8f\n", guideCurves(density, :, 27));
fprintf(fileID, "END\t%1f\n", 0); 

for (i = 1:guideCurveCount)
    fprintf(fileID, "GuideCurve%1.0f\n", i);
    fprintf(fileID, "START\t%1f\n", 0); 
    for(j = 1:density)    
        fprintf(fileID, "%9.8f\t%9.8f\t%9.8f\n", guideCurves(j, :, i));
    end
    fprintf(fileID, "END\t%1f\n", 0); 
end
fprintf(fileID, "HubRadius %3.5f\n", rh);
fprintf(fileID, "Blades %1.f\n", Z);
fclose(fileID);

