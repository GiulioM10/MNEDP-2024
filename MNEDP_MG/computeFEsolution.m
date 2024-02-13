function [U, u] = computeFEsolution(geom, A, b, gD_f, iterative)

arguments
    geom
    A
    b
    gD_f
    iterative = true
end
%UNTITLED3 Solves the linear system and applies Dirichlet conditions
%   Detailed explanation goes here

% Solve the linear system
perm = symrcm(A);
if iterative
    try R = ichol(A(perm, perm), struct('michol','on'));
        disp("Using pcg...")
        [z,flag,relres] = pcg(A(perm, perm), b(perm), 1e-8, 800, R, R');
    catch
        disp("Using gmres...")
        [L, U] = ilu(A(perm, perm));
        [z, flag, relres] = gmres(A(perm, perm), b(perm), 100, 1e-8, 20, L, U);
    end
    if flag
        warning("Converegnce was not reached")
        disp(relres)
    end
else
    disp("Using \...")
    z = A(perm, perm)\b(perm);
end
u(perm) = z;

% Allocate memory for the complete solution vector
U = zeros(geom.Nobj.N_node, 1);

for i = 1:geom.Nobj.N_node % For every node in the triangulation ...
    ii = geom.piv.piv(i);
    if ii > 0 % ... if it is a degree of freedom, assign the numerical value
        U(i) = u(geom.piv.piv(i));
    else % ... otherwise, assign the value given by the Dirichlet condition
        U(i) = gD_f(geom.obj.P(i, 1), geom.obj.P(i, 2));
    end
end
end