function [err_L2, err_H1, err_Linf] = computeError(k, geom, U, Utrue, gradUtrue)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[Phi_hat_f, Dx_Phi_hat_f, Dy_Phi_hat_f] = getBasisFunctions(k);

[xhat, yhat, omega, numNodes] = getNodesWeights();
Phi_hat = Phi_hat_f(xhat, yhat);
Dx_Phi_hat = Dx_Phi_hat_f(xhat, yhat);
Dy_Phi_hat = Dy_Phi_hat_f(xhat, yhat);

dx = zeros(2, 1);
dy = zeros(2, 1);
err_L2 = 0.0;
err_H1 = 0.0;

for e = 1:geom.Nobj.N_ele % For each element...
    % dx_1 = x_3 - x_2
    dx(1) = geom.obj.P(geom.obj.T(e, 3), 1)...
        - geom.obj.P(geom.obj.T(e, 2), 1);
    % dx_2 = x_1 - x_3
    dx(2) = geom.obj.P(geom.obj.T(e, 1), 1)...
        - geom.obj.P(geom.obj.T(e, 3), 1);
    % dy_1 = y_2 - y_3
    dy(1) = geom.obj.P(geom.obj.T(e, 2), 2)...
        - geom.obj.P(geom.obj.T(e, 3), 2);
    % dy_2 = y_3 - y_1
    dy(2) = geom.obj.P(geom.obj.T(e, 3), 2)...
        - geom.obj.P(geom.obj.T(e, 1), 2);
    
    % F_e (xhat) = B*xhat + b
    B = [dx(2), -dx(1);
         -dy(2), dy(1)];
    b = geom.obj.P(geom.obj.T(e, 3), :)';

    Btinv = 1/(dx(2)*dy(1) - dx(1)*dy(2)) * ...
        [dy(1), dy(2);
         dx(1), dx(2)];

    % Retrieve the area of the element
    Area = geom.support.TInfo(e).Area;
    

    for q = 1:numNodes % For each node...
        x = b + B * [xhat(q); yhat(q)]; % Map the node to the element
        % Update the L2 error
        err_L2 = err_L2 + omega(q) * 2 * Area * ...
            abs(Utrue(x(1), x(2)) - Phi_hat(q, :)*U(geom.obj.T(e, :)))^2;
        % Update the H1 error
        err_H1 = err_H1 + omega(q) * 2 * Area * ...
            (gradUtrue(x(1), x(2)) - Btinv * [Dx_Phi_hat(q, :); Dy_Phi_hat(q, :)]*U(geom.obj.T(e, :)))'*...
            (gradUtrue(x(1), x(2)) - Btinv * [Dx_Phi_hat(q, :); Dy_Phi_hat(q, :)]*U(geom.obj.T(e, :)));
    end
end
% Compute the square root
err_L2 = sqrt(err_L2);
err_H1 = sqrt(err_H1);

u_true_vec = zeros(geom.Nobj.N_node, 1);
for i = 1:geom.Nobj.N_node
    u_true_vec(i) = Utrue(geom.obj.P(i, 1), geom.obj.P(i, 2));
end

err_Linf = norm(U-u_true_vec, "inf");

end