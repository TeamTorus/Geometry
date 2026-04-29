function [cartX, cartY, cartZ] = Wrapping3D(airfoilBase, workingRadius, scalingFactor, attack, skew, rake)
    %airfoilBase is a Nx2 list of (x, y) coordinates of the airfoil
    %workingRadius is a scalar of the current radius
    %scalingFactor is a scalar of the scale being applied to the airfoil
    %attack, skew, and rake are scalars

    airfoilBase = airfoilBase .* scalingFactor; %scale the airfoil
    airfoilBase = airfoilBase * [cosd(attack), -sind(attack); sind(attack), cosd(attack)]; %rotate the airfoil by the AoA
    z = airfoilBase(:, 2) - rake; %offset all of the y points by the rake and represent them as the z coordinate
    theta = (airfoilBase(:, 1) ./ workingRadius) - skew; %find the wrapping theta and apply skew
    [cartX, cartY, cartZ] = pol2cart(theta, workingRadius, z); %transform from polar to cartesian
end