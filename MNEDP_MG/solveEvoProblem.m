function [u, U, condA] = solveEvoProblem(K, geom, t, Dt, rho, mu, beta, sigma, f, gD_f, Dt_gD_f, gN_f, u0_f, SUPG, MassLumping, iterative, forcePcg, forceGmres)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
arguments
    K 
    geom struct
    t
    Dt
    rho function_handle
    mu function_handle
    beta function_handle
    sigma function_handle
    f function_handle
    gD_f function_handle
    Dt_gD_f function_handle
    gN_f function_handle
    u0_f function_handle
    SUPG logical = false
    MassLumping logical = false
    iterative logical = true
    forcePcg logical = false
    forceGmres logical = false
end
%%% EXTRACT USEFUL DATA FROM geom %%%
Ndof = max(geom.piv.piv); % # of degerees of freedom
ND = size(geom.piv.Di(:, 1), 1); % # of nodes with Dirichlet conds.
NN = size(geom.piv.Ne, 1); % # of edges with Neuman conds.
NVert = (K + 2)*(K + 1)/2; % # of nodes on each element
Nt = size(t, 2); % # of time instants


%%% DATA STRUCTURES FOR EFFICIENTLY BUILDING SPARSE MATRICES %%%
% Stifness Matrix %
rows = zeros(10 * Ndof, 1); % Vector containing row indices
cols = zeros(10 * Ndof, 1); % Vector containing column indices
values = zeros(10 * Ndof, 1); % Vector containing values of the matrix
count = 0; % Position of the last written element in the vectors

% Dirichlet Matrix %
rowsD = zeros(5 * Ndof, 1); % Vector containing row indices
colsD = zeros(5 * Ndof, 1); % Vector containing column indices
valuesD = zeros(5 * Ndof, 1); % Vector containing values of the matrix
countD = 0; % Position of the last written element in the vectors

% Mass Matrix %
rowsM = zeros(10 * Ndof, 1); % Vector containing row indices
colsM = zeros(10 * Ndof, 1); % Vector containing column indices
valuesM = zeros(10 * Ndof, 1); % Vector containing values of the matrix
countM = 0; % Position of the last written element in the vectors

% Mass Dirichlet Matrix %
rowsMD = zeros(5 * Ndof, 1); % Vector containing row indices
colsMD = zeros(5 * Ndof, 1); % Vector containing column indices
valuesMD = zeros(5 * Ndof, 1); % Vector containing values of the matrix
countMD = 0; % Position of the last written element in the vectors

b = zeros(Ndof, 1); % Right hand-term relative to time t_{n+1}
gD = zeros(ND, 1); % "Rilevamento" vector
Dt_gD = zeros(ND, 1); % Time derivative of "rilevamento"
b_prev = zeros(Ndof, 1);  % Right hand-term relative to time t_n
u_prev = zeros(Ndof, 1); % u vector relative to time t_n

%%% PREPARE VARIABLES FOR COMPUTING INTEGRALS %%%
[Phi_hat_f, Dx_Phi_hat_f, Dy_Phi_hat_f] = getBasisFunctions(K);

% Get the nodes and the weights for integral computation
[xhat, yhat, omega, numNodes] = getNodesWeights();

% Evaluate the basis function in the nodes
Phi_hat = Phi_hat_f(xhat, yhat);
Dx_Phi_hat = Dx_Phi_hat_f(xhat, yhat);
Dy_Phi_hat = Dy_Phi_hat_f(xhat, yhat);

if SUPG % If the user requires a diffusion-convection stabilization...
    %... compute the following quantities

    %%% CONSTANT OF INVERSE DISEQUALITY %%%
    mk = 1/24 * (K == 2) + 1/3 * (K == 1);

    %%% SECOND DERIVATIVES OF BASIS FUNCTIONS %%%
        D2Phi = [
            4, 0, 4, -8, 0,  0;
            0, 0, 4, -4, 4, -4;
            0, 4, 4,  0, 0, -8;
            ] * (K - 1);
end %if SUPG

%%% TIME INDEPENDENT PART + RIGHT HAND TERM AT TIME t_0 %%%
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
    
    % x = B * x_hat + a3
    B = [dx(2), -dx(1);
         -dy(2), dy(1)];
    a3 = geom.obj.P(geom.obj.T(e, 3), :)';
    
    % Manually compute the inverse of B^t
    Btinv = 1/(dx(2)*dy(1) - dx(1)*dy(2)) * ...
        [dy(1), dy(2);
         dx(1), dx(2)];

    % Retrieve the area of the element
    Area = geom.support.TInfo(e).Area;
    
    if SUPG
        % Estimate the diameter of the element
        hE = sqrt(Area);
    
        % Retrieve the center of mass of the element
        G = geom.support.TInfo(e).CG;
    
        % Compute Peclè number of the element...
        normBeta = norm(beta(G(1), G(2)));
        epsilon = mu(G(1), G(2));
        Peh = mk * normBeta * hE / (2 * epsilon);
    
        % ... and virtual viscosity accordingly
        if Peh <= 1
            tauE = mk * hE^2 / (4 * epsilon);
        else
            tauE = hE / (2 * normBeta);
        end
    end
    
    for j = 1:NVert % For each node on the element...
        jj = geom.piv.piv(geom.obj.T(e, j)); % <=> Pivot(Global(j))

        if jj > 0 % If it is a degree of freedom...

            if MassLumping % If the user requires a diffusion-reaction stabilization...
                reactionLumpTerm = 0.0; % Prepare lumped term
            end

            for k = 1:NVert % For each node...
                kk = geom.piv.piv(geom.obj.T(e, k)); % <=> Pivot(Global(k))

                % Prepare variables for integral computation
                diffusionTerm = 0.0;
                convectionTerm = 0.0;
                reactionTerm = 0.0;
                stabilizationTerm = 0.0;

                massTerm = 0.0;
                
                % Compute integrals
                for q = 1:numNodes % For each integration node...
                    % x = F_e(x_hat)
                    x = a3 + B * [xhat(q); yhat(q)];
                    
                    % d_jk = int_E(mu * dot(grad(phi_k), grad(phi_j)))
                    diffusionTerm = diffusionTerm + omega(q) * 2 * Area * ...
                        mu(x(1), x(2)) * ...
                        (Btinv * [Dx_Phi_hat(q, k); Dy_Phi_hat(q, k)])' * ...
                        (Btinv * [Dx_Phi_hat(q, j); Dy_Phi_hat(q, j)]);
                    
                    % c_jk = int_E(dot(beta, grad(phi_k)) * phi_j)
                    convectionTerm = convectionTerm + omega(q) * 2 * Area * ...
                        beta(x(1), x(2))' * ...
                        Btinv * [Dx_Phi_hat(q, k); Dy_Phi_hat(q, k)] * Phi_hat(q, j);

                    % r_jk = int_E(sigma * phi_k * phi_j)
                    reactionTerm = reactionTerm + omega(q) * 2 * Area * ...
                        sigma(x(1), x(2)) * Phi_hat(q, k) * Phi_hat(q, j);
                    
                    % m_jk = int_E(phi_k*phi_j)
                    massTerm = massTerm + omega(q) * 2 * Area * ...
                        rho(x(1), x(2)) * Phi_hat(q, k) * Phi_hat(q, j);

                    if SUPG % If the user requires a diffusion-convection stabilization...
                        % s_jk = int_E(tau_E * (laplacian(phi_k) + dot(beta, grad(phi_k)))*dot(beta, grad(phi_j))
                        stabilizationTerm = stabilizationTerm + omega(q) * 2 * Area * ...
                            tauE*(beta(x(1), x(2))' * ...
                            Btinv * [Dx_Phi_hat(q, k); Dy_Phi_hat(q, k)] ...
                            - mu(x(1), x(2)) * [dy(1)^2 + dx(1)^2, 2*(dy(1)*dy(2) + dx(1)*dx(2)), dy(2)^2 + dx(2)^2] * D2Phi(:, k)) * ...
                            beta(x(1), x(2))' * Btinv * [Dx_Phi_hat(q, j); Dy_Phi_hat(q, j)];
                    end
                end % q
                if kk > 0 % If it is a degree of freedom...
                    % ... update the stifness matrix ...
                    count = count + 1; % Write in a new cell
                    rows(count) = jj; % Save row index
                    cols(count) = kk; % Save column index
                    values(count) = diffusionTerm + convectionTerm + stabilizationTerm + reactionTerm * (~MassLumping); % Save value
                    
                    if MassLumping
                        reactionLumpTerm = reactionLumpTerm + reactionTerm; % Add the contribution to lumped term
%                         massLumpTerm = massLumpTerm + massTerm;
                    end

                    % ... update the mass matrix ...
                    countM = countM + 1; % Write in a new cell
                    rowsM(countM) = jj; % Save row index
                    colsM(countM) = kk; % Save column index
                    valuesM(countM) = massTerm; % Save value

                else % ... otherwise ...
                    % ... update the Dirichlet matrix
                    countD = countD + 1; % Write in a new cell
                    rowsD(countD) = jj; % Save row index
                    colsD(countD) = -kk; % Save column index
                    valuesD(countD) = diffusionTerm + convectionTerm + stabilizationTerm + reactionTerm; % Save value

                    % ... update the Mass Dirichlet matrix
                    countMD = countMD + 1; % Write in a new cell
                    rowsMD(countMD) = jj; % Save row index
                    colsMD(countMD) = -kk; % Save column index
                    valuesMD(countMD) = massTerm; % Save value
                end %kk
            end % k

            if MassLumping
                % Write contribution to the appropriate cell
                count = count + 1; % Write in a new cell
                rows(count) = jj; % Save row index
                cols(count) = jj; % Save column index
                values(count) = reactionLumpTerm; % Save value
            end
            % Prepare variables for integral computation
            fTerm = 0.0;
            fstabilizationTerm = 0.0;
            
            % Compute integral
            for q = 1:numNodes % For each integration node...
                % x = F_e(x_hat)
                x = a3 + B*[xhat(q); yhat(q)]; 
                
                % f_j = int_E(f*phi_j)
                fTerm = fTerm + omega(q) * 2 * Area * f(x(1), x(2), t(1)) * Phi_hat(q, j);
                
                if SUPG
                    % fs_j = int_E(f*dot(beta, grad(phi_j)))
                    fstabilizationTerm = fstabilizationTerm + omega(q) * 2 * Area * tauE * ...
                        f(x(1), x(2), t(1)) * beta(x(1), x(2))' * Btinv * [Dx_Phi_hat(q, j); Dy_Phi_hat(q, j)];
                end
            end % q
            % Update the known term
            b_prev(jj) = b_prev(jj) + fTerm + fstabilizationTerm;
        end % jj
    end % j
end % e

% Trim vectors and build Stifness Matrix
rows = rows(1:count);
cols = cols(1:count);
values = values(1:count);
A = sparse(rows, cols, values, Ndof, Ndof);

condA = condest(A);

% Trim vectors and build Dirichlet Matrix
rowsD = rowsD(1:countD);
colsD = colsD(1:countD);
valuesD = valuesD(1:countD);
AD = sparse(rowsD, colsD, valuesD, Ndof, ND);

% Trim vectors and build Mass Matrix
rowsM = rowsM(1:countM);
colsM = colsM(1:countM);
valuesM = valuesM(1:countM);
M = sparse(rowsM, colsM, valuesM, Ndof, Ndof);

% Trim vectors and build Mass Dirichlet Matrix
rowsMD = rowsMD(1:countMD);
colsMD = colsMD(1:countMD);
valuesMD = valuesMD(1:countMD);
MD = sparse(rowsMD, colsMD, valuesMD, Ndof, ND);

% Crank-Nicolson's method matrices
Z = M + (Dt/2)*A;
W = M - (Dt/2)*A;

solver = 0;
perm = symrcm(Z);
if iterative
    if forcePcg
        R = ichol(Z(perm, perm), struct('michol','on'));
        disp("Using pcg...")
        solver = 1;
    elseif forceGmres
        disp("Using gmres...")
        [L, U] = ilu(Z(perm, perm));
        solver = 2;
    else
    try R = ichol(Z(perm, perm), struct('michol','on'));
        disp("Using pcg...")
        solver = 1;
    catch
        disp("Using gmres...")
        [L, U] = ilu(Z(perm, perm));
        solver = 2;
    end
    end
else
    disp("Using \...")
end

% Neumann conditions at time t_0
switch K
    case 1 % If using P_1...
        NMatrix = [1/3, 1/6;
                    1/6, 1/3]; % Useful matrix to speed up calculation

        for e = 1:NN % For every edge on the Neumann border...
            % Extract the beginning point the end point and compute the length of
            % the edge
            xb = geom.obj.P(geom.obj.E(geom.piv.Ne(e, 1), 1), :); 
            xe = geom.obj.P(geom.obj.E(geom.piv.Ne(e, 1), 2), :);
            len_e = norm(xe - xb, 2);
            
            % Evaluate g_N in those points
            gN = [gN_f(xb(1), xb(2), t(1)); gN_f(xe(1), xe(2), t(1))];

            for i = 1:2 % For each extreme point
                j = geom.piv.piv(geom.obj.E(geom.piv.Ne(e, 1), i));
                if j > 0 % If it is a degree of freedom
                    b_prev(j) = b_prev(j) + NMatrix(i, :) * gN * len_e; % Update the term
                end
            end
        end
    case 2 % If using P_2...
        NMatrix = [2/15, 1/15, -1/30;
                   1/15, 8/15,  1/15;
                  -1/30, 1/15,  2/15]; % Useful matrix to speed up calculation

        for e = 1:NN % For every edge on the Neumann border...
            % Extract the beginning, middle and end points and compute the length of
            % the edge
            xb = geom.obj.P(geom.obj.E(geom.piv.Ne(e, 1), 1), :); 
            xe = geom.obj.P(geom.obj.E(geom.piv.Ne(e, 1), 2), :);
            xm = geom.obj.P(geom.obj.E(geom.piv.Ne(e, 1), 5), :);
            len_e = norm(xe - xb, 2);
            
            % Evaluate gN in those points
            gN = [gN_f(xb(1), xb(2), t(1)); gN_f(xm(1), xm(2), t(1)); gN_f(xe(1), xe(2), t(1))];
            idx = [1, 5, 2]; % Index in wich are stored info relative to each point

            for i = 1:3 % For each node...
                j = geom.piv.piv(geom.obj.E(geom.piv.Ne(e, 1), idx(i)));
                if j > 0 % If it is a degree of freedom...
                    b_prev(j) = b_prev(j) + NMatrix(i, :) * gN * len_e; % Update the term
                end
            end
        end
    otherwise % Other methods are not implemented
        error("You shouldn't be here")
end


% Dirichlet conditions at time t_0
% Compute the vector containing all the vaues of the "rilevamento"
for i = 1:ND
    ii = geom.piv.Di(i, 1);
    gD(i) = gD_f(geom.obj.P(ii, 1), geom.obj.P(ii, 2), t(1));

    Dt_gD(i) = Dt_gD_f(geom.obj.P(ii, 1), geom.obj.P(ii, 2), t(1));
    % Derivatives are approximated using finite differences
end
% Update the right at time t_0 hand term accordingly
b_prev = b_prev - AD * gD - MD * Dt_gD;

% Initialize vector containing values of initial condition
for i = 1:geom.Nobj.N_node
    ii = geom.piv.piv(i);
    if ii > 0
        u_prev(ii) = u0_f(geom.obj.P(i, 1), geom.obj.P(i, 2));
    end
end

%%% CRANK-NICKOLSON METHOD %%%
for n = 2:Nt % for t_1, t_2, ... t_N
% Assemble right hand term
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
    
        % x = B * x_hat + a3
        B = [dx(2), -dx(1);
            -dy(2), dy(1)];
        a3 = geom.obj.P(geom.obj.T(e, 3), :)';
    
        % Manually compute the inverse of B^t
        Btinv = 1/(dx(2)*dy(1) - dx(1)*dy(2)) * ...
            [dy(1), dy(2);
            dx(1), dx(2)];

        % Retrieve the area of the element
        Area = geom.support.TInfo(e).Area;
    
        if SUPG
            % Estimate the diameter of the element
            hE = sqrt(Area);
    
            % Retrieve the center of mass of the element
            G = geom.support.TInfo(e).CG;
    
            % Compute Peclè number of the element...
            normBeta = norm(beta(G(1), G(2)));
            epsilon = mu(G(1), G(2));
            Peh = mk * normBeta * hE / (2 * epsilon);
    
            % ... and virtual viscosity accordingly
            if Peh <= 1
                tauE = mk * hE^2 / (4 * epsilon);
            else
                tauE = hE / (2 * normBeta);
            end
        end

        for j = 1:NVert % For each node on the element...
            jj = geom.piv.piv(geom.obj.T(e, j)); % <=> Pivot(Global(j))

            if jj > 0 % If it is a degree of freedom...
                % Prepare variables for integral computation
                fTerm = 0.0;
                fstabilizationTerm = 0.0;
            
                % Compute integral
                for q = 1:numNodes % For each integration node...
                    % x = F_e(x_hat)
                    x = a3 + B*[xhat(q); yhat(q)]; 
                
                    % f_j = int_E(f*phi_j)
                    fTerm = fTerm + omega(q) * 2 * Area * f(x(1), x(2), t(n)) * Phi_hat(q, j);
                
                    if SUPG
                        % fs_j = int_E(f*dot(beta, grad(phi_j)))
                        fstabilizationTerm = fstabilizationTerm + omega(q) * 2 * Area * tauE * ...
                            f(x(1), x(2), t(n)) * beta(x(1), x(2))' * Btinv * [Dx_Phi_hat(q, j); Dy_Phi_hat(q, j)];
                    end
                end % q
                % Update the known term
                b(jj) = b(jj) + fTerm + fstabilizationTerm;
            end % jj
        end % j
    end %e

    % Dirichlet conditions
    % Compute the vector containing all the vaues of the "rilevamento"
    for i = 1:ND
        ii = geom.piv.Di(i, 1);
        gD(i) = gD_f(geom.obj.P(ii, 1), geom.obj.P(ii, 2), t(n));
        Dt_gD(i) = Dt_gD_f(geom.obj.P(ii, 1), geom.obj.P(ii, 2), t(n));
    end
    b = b - AD * gD - MD * Dt_gD;

    % Neumann condtions
    switch K
        case 1
            NMatrix = [1/3, 1/6;
                        1/6, 1/3]; % Useful matrix to speed up calculation

            for e = 1:NN % For every edge on the Neumann border.....
                % Extract the beginning point the end point and compute the length of
                % the edge
                xb = geom.obj.P(geom.obj.E(geom.piv.Ne(e, 1), 1), :); 
                xe = geom.obj.P(geom.obj.E(geom.piv.Ne(e, 1), 2), :);
                len_e = norm(xe - xb, 2);
                gN = [gN_f(xb(1), xb(2), t(n)); gN_f(xe(1), xe(2), t(n))];
                for i = 1:2 % For each extreme point
                    j = geom.piv.piv(geom.obj.E(geom.piv.Ne(e, 1), i));
                    if j > 0 % If it is a degree of freedom
                        b(j) = b(j) + NMatrix(i, :) * gN * len_e; % Update the term
                    end
                end
            end
        case 2
            NMatrix = [2/15, 1/15, -1/30;
                        1/15, 8/15,  1/15;
                        -1/30, 1/15,  2/15]; % Useful matrix to speed up calculation

            for e = 1:NN
                xb = geom.obj.P(geom.obj.E(geom.piv.Ne(e, 1), 1), :); 
                xe = geom.obj.P(geom.obj.E(geom.piv.Ne(e, 1), 2), :);
                xm = geom.obj.P(geom.obj.E(geom.piv.Ne(e, 1), 5), :);
                gN = [gN_f(xb(1), xb(2), t(n)); gN_f(xm(1), xm(2), t(n)); gN_f(xe(1), xe(2), t(n))];
                len_e = norm(xe - xb, 2);
                idx = [1, 5, 2];
                for i = 1:3
                    j = geom.piv.piv(geom.obj.E(geom.piv.Ne(e, 1), idx(i)));
                    if j > 0
                        b(j) = b(j) + NMatrix(i, :) * gN * len_e;
                    end
                end
            end
        otherwise
            error("You shouldn't be here")
    end
    
    % Generate permutation to prevent eccessive fill-in

    if solver == 0
        u_temp = Z(perm, perm)\(W(perm, perm)*u_prev(perm) + (Dt/2) *(b_prev(perm) + b(perm)));
    else
        if solver == 1
            [u_temp,flag,relres] = pcg(Z(perm, perm), W(perm, perm)*u_prev(perm) + (Dt/2) *(b_prev(perm) + b(perm)), 1e-8, 1500, R, R');
        elseif solver == 2
            [u_temp, flag, relres] = gmres(Z(perm, perm), W(perm, perm)*u_prev(perm) + (Dt/2) *(b_prev(perm) + b(perm)), 100, 1e-8, 20, L, U);
        else
            error("Something is worng")
        end
        if flag
            warning("Converegnce was not reached")
            disp(relres)
        end
    end

    u_prev(perm) = u_temp;
    b_prev = b;
    b = zeros(Ndof, 1); 
end %n

% Allocate memory for the complete solution vector
U = zeros(geom.Nobj.N_node, 1);

u = u_prev;
for i = 1:geom.Nobj.N_node % For every node in the triangulation ...
    ii = geom.piv.piv(i);
    if ii > 0 % ... if it is a degree of freedom, assign the numerical value
        U(i) = u(ii);
    else % ... otherwise, assign the value given by the Dirichlet condition
        U(i) = gD_f(geom.obj.P(i, 1), geom.obj.P(i, 2), t(end));
    end
end
end