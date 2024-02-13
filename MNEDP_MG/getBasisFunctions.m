function [Phi_hat_f, Dx_Phi_hat_f, Dy_Phi_hat_f] = getBasisFunctions(k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N1 = @(x, y) x;
N2 = @(x, y) y;
N3 = @(x, y) 1 - x - y;

dx_N1 = @(x, y) ones(size(x));
dx_N2 = @(x, y) zeros(size(x));
dx_N3 = @(x, y) -ones(size(x));

dy_N1 = @(x, y) zeros(size(x));
dy_N2 = @(x, y) ones(size(x));
dy_N3 = @(x, y) -ones(size(x));

switch k
    case 1
        Phi_hat_f = @(x, y) ...
            [N1(x, y), N2(x, y), N3(x, y)];
        Dx_Phi_hat_f = @(x, y) ...
            [dx_N1(x, y), dx_N2(x, y), dx_N3(x, y)];
        Dy_Phi_hat_f = @(x, y) ...
            [dy_N1(x, y), dy_N2(x, y), dy_N3(x, y)];

    case 2
        Phi_hat_f = @(x, y) ...
            [2 * N1(x, y).*(N1(x, y) - .5), 2 * N2(x, y).*(N2(x, y) - .5), 2 * N3(x, y).*(N3(x, y) - .5), ...
            4 * N1(x, y).*N3(x, y), 4 * N1(x, y).*N2(x, y), 4 * N2(x, y).*N3(x, y)];

        Dx_Phi_hat_f = @(x, y) ...
            [2 * dx_N1(x, y).*(N1(x, y) - .5) + 2 * N1(x, y).*dx_N1(x, y), 2 * dx_N2(x, y).*(N2(x, y) - .5) + 2 * N2(x, y).*dx_N2(x, y), 2 * dx_N3(x, y).*(N3(x, y) - .5) + 2 * N3(x, y).*dx_N3(x, y), ...
            4 * dx_N1(x, y).*N3(x, y) + 4 * N1(x, y).*dx_N3(x, y), 4 * dx_N1(x, y).*N2(x, y) + 4 * N1(x, y).*dx_N2(x, y), 4 * dx_N2(x, y).*N3(x, y) + 4 * N2(x, y).*dx_N3(x, y)];

        Dy_Phi_hat_f = @(x, y) ...
            [2 * dy_N1(x, y).*(N1(x, y) - .5) + 2 * N1(x, y).*dy_N1(x, y), 2 * dy_N2(x, y).*(N2(x, y) - .5) + 2 * N2(x, y).*dy_N2(x, y), 2 * dy_N3(x, y).*(N3(x, y) - .5) + 2 * N3(x, y).*dy_N3(x, y), ...
            4 * dy_N1(x, y).*N3(x, y) + 4 * N1(x, y).*dy_N3(x, y), 4 * dy_N1(x, y).*N2(x, y) + 4 * N1(x, y).*dy_N2(x, y), 4 * dy_N2(x, y).*N3(x, y) + 4 * N2(x, y).*dy_N3(x, y)];

    otherwise
        error("Finite element method based on P" + int2str(k) + " not implemented yet.");
end
end