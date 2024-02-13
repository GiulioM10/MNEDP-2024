clear all
close all
clc

if(~exist('mybbtr30.m'))
     addpath('../bbtr30')
     addpath('../MNEDP_MG')
     disp('../bbtr30 added to the path')
end

%% DEFINE THE DOMAIN
dVertices = [ 0 0
              1 0
              1 1
              0 1];
dBoundary = 1:4;
bcBoundary = [9 11 2 13];
bcVertices = [1 3 5 7];
bcValues = [3.0, 2.0];
checkArea = "Y";
checkAngle = "N";
areaValue = 0.001;
angleValue = 30;

t0 = 0;
T = 1;

K = 2;

%%
utrue_t = @(x, y, t) t^2 + sin(2*pi*(x + y));

dt_utrue_t = @(x, y, t) 2*t;

gradientutrue_t = @(x, y, t) [2*pi*cos(2*pi*(x + y))
                              2*pi*cos(2*pi*(x + y))];

utrue = @(x, y) 1 + sin(2*pi*(x + y));

gradientutrue = @(x, y) [2*pi*cos(2*pi*(x + y))
                         2*pi*cos(2*pi*(x + y))];
 
laplacianutrue = @(x, y, t) 2*x*pi*cos(2*pi*(x + y))*exp(x*y) - 8*pi^2*exp(x*y)*sin(2*pi*(x + y)) + 2*y*pi*cos(2*pi*(x + y))*exp(x*y);

rho = @(x, y) exp(-x -y);

mu = @(x, y) exp(x*y);

beta = @(x, y) [-x; y];

sigma = @(x, y) -x-y;

f = @(x, y, t) rho(x, y)*dt_utrue_t(x, y, t) -laplacianutrue(x, y, t) + beta(x, y)'*gradientutrue_t(x,y,t) + sigma(x, y)*utrue_t(x, y, t);

gD_f = @(x, y, t) t^2 + sin(2*pi*(x + y));

Dt_gD_f = @(x, y, t) 2*t;

gN_f = @(x, y, t) mu(x, y) .* [0, 1] * gradientutrue_t(x, y, t);

u0_f = @(x, y) sin(2*pi*(x + y));

%%
draw = true;
areaValue = .01*[0.08, 0.04, 0.01];
errorVector = zeros(6, length(areaValue));
[t, dt] = discretizeTime(t0, T, 0.05);

for exp = 1:length(areaValue)
    geom = defineTriangulation(K, dVertices, dBoundary, bcBoundary, ...
        bcVertices, bcValues, checkArea, checkAngle, areaValue(exp), angleValue, false);

    [u, U, condA] = solveEvoProblem(K, geom, t, dt, rho, mu, beta, sigma, f, gD_f, Dt_gD_f, gN_f, u0_f);

    [err_L2, err_H1, err_Linf] = computeError(K, geom, U, utrue, gradientutrue);

    errorVector(1, exp) = sqrt(max([geom.support.TInfo(:).Area]));
    errorVector(2, exp) = err_L2;
    errorVector(3, exp) = err_H1;
    errorVector(4, exp) = err_Linf;
    errorVector(5, exp) = sqrt(min([geom.support.TInfo(:).Area]));
    errorVector(6, exp) = condA;

    if exp == 1 && draw
        x = linspace(0, 1);
        y = linspace(0, 1);
        [X, Y] = meshgrid(x, y);
        Z = utrue(X, Y);

        printComparison(K, geom, U, X, Y, Z)
    end
end

% approxDt =  [0.1, 0.08, 0.06];
% errorVectorT = zeros(3, length(approxDt));
% geom = defineTriangulation(K, dVertices, dBoundary, bcBoundary, ...
%     bcVertices, bcValues, checkArea, checkAngle, areaValue(1), angleValue, false);
% 
% for exp = 1:length(approxDt)
%     [t, dt] = discretizeTime(t0, T, approxDt(exp));
% 
%     [u, U, ~] = solveEvoProblem(K, geom, t, dt, rho, mu, beta, sigma, f, gD_f, Dt_gD_f, gN_f, u0_f);
% 
%     [err_L2, err_H1, ~] = computeError(K, geom, U, utrue, gradientutrue);
% 
%     errorVectorT(1, exp) = dt;
%     errorVectorT(2, exp) = err_L2;
%     errorVectorT(3, exp) = err_H1;
% end



figure
tiledlayout(2, 3)
nexttile
semilogy(errorVector(1, :), errorVector(2, :), "b-diamond")
title("Convergenza in L^2 (h)")
xlabel("h_{max}")
ylabel("log(E_0)")
grid on
nexttile
semilogy(errorVector(1, :), errorVector(3, :), "m-diamond")
title("Convergenza in H^1")
xlabel("h_{max}")
ylabel("log(E_1)")
grid on
nexttile
plot(errorVector(1, :), errorVector(6, :), "k-diamond")
title("Condizionamento di A")
xlabel("h_{min}")
ylabel("cond_2(A)")
grid on
nexttile([1,3])
semilogy(errorVector(1, :), errorVector(2, :), "b-diamond", errorVector(1, :), errorVector(4, :), "r-diamond")
title("Ordini di convergenza")
xlabel("h_{max}")
legend("log(E_0)", "log(E_\infty)", 'Location', 'northwest')
grid on

% figure
% tiledlayout(1, 2)
% nexttile
% semilogy(errorVectorT(1, :), errorVectorT(2, :), "b-diamond")
% title("Convergenza in L^2 (\Delta t)")
% xlabel("\Delta t")
% grid on
% nexttile
% semilogy(errorVectorT(1, :), errorVectorT(3, :), "m-diamond")
% title("Convergenza in H^1 (\Delta t)")
% xlabel("\Delta t")
% ylabel("log(E_1)")
% grid on

clc
l2 = polyfit(log(errorVector(1, :)), log(errorVector(2, :)), 1);
l2 = l2(1);
h1 = polyfit(log(errorVector(1, :)), log(errorVector(3, :)), 1);
h1 = h1(1);
linf = polyfit(log(errorVector(1, :)), log(errorVector(4, :)), 1);
linf = linf(1);
cond = polyfit(log(errorVector(5, :)), log(errorVector(6, :)), 1);
cond = cond(1);
% l2t = polyfit(log(errorVectorT(1, :)), log(errorVectorT(2, :)), 1);
% l2t = l2t(1);
% h1t = polyfit(log(errorVectorT(1, :)), log(errorVectorT(3, :)), 1);
% h1t = h1t(1);
fprintf("------------------------------- \n")
fprintf("Ordine di convergenza in L2 = %d \n", round(l2))
fprintf("Ordine di convergenza in H1 = %d \n", round(h1))
fprintf("Ordine di convergenza in LInf = %d \n", round(linf))
fprintf("Ordine di crescita di cond(A) = %d \n", round(cond))
% fprintf("Ordine di convergenza in L2 (dt) = %d \n", round(l2t))
% fprintf("Ordine di convergenza in H1 (dt) = %d \n", round(h1t))
fprintf("------------------------------- \n")

