function [x, y] = NACA66TMBModBuilder(chordLength, thicknessMax, camberMax)
    XTCD = importdata ("NACA66TMBModXTCD.txt");
    xNorm = XTCD(:, 1) * chordLength;
    yNorm = XTCD(:, 3) * camberMax;
    tNorm = XTCD(:, 2) * thicknessMax;
    dNorm = XTCD(:, 4) * camberMax;
    x = xNorm - tNorm .* sin(atan(dNorm));
    y = yNorm + tNorm .* cos(atan(dNorm));
end