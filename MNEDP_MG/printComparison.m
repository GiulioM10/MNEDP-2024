function [] = printComparison(k, geom,U, X, Y, Z, drawlines)
%UNTITLED4 Print a side by side comparison of the numerical solution with
%the real one
%   Detailed explanation goes here

arguments
    k
    geom
    U
    X
    Y
    Z
    drawlines logical = true
end

figure % New graphical window
tiledlayout(1, 2) % Multiplot

nexttile
switch k
    case 1
        T = geom.obj.T;
    case 2
        T = plottingTraingulation(geom.obj.T);
    otherwise
        error("Finite element method based on P" + int2str(k) + " not implemented yet.");
end
if drawlines
    trisurf(T, geom.obj.P(:, 1), geom.obj.P(:, 2), U); % plot of u_h
else
    h = trisurf(T, geom.obj.P(:, 1), geom.obj.P(:, 2), U); % plot of u_h
    set(h, 'edgecolor', 'none')
end

nexttile
surf(X, Y, Z) % plot of u

figure()
if drawlines
    trisurf(T, geom.obj.P(:, 1), geom.obj.P(:, 2), U); % plot of u_h
else
    h = trisurf(T, geom.obj.P(:, 1), geom.obj.P(:, 2), U); % plot of u_h
    set(h, 'edgecolor', 'none')
end
figure()
surf(X, Y, Z) % plot of u
end