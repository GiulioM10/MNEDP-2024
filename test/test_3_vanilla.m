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

% approxDt =  [0.08, 0.06, 0.04];
% t0 = 0;
% T = 1;

K = 2;

%%

utrue = @(x, y) tanh(40*(x - y))+1;

gradientutrue = @(x, y) [40 - 40*tanh(40*x - 40*y)^2;
                     40*tanh(40*x - 40*y)^2 - 40];
 
laplacianutrue = @(x,y) - exp(x*y)*(2000*x*(x - 1)*(x - 1/2)*(y - 1) + 2000*x*(x - 1)*(x - 1/2)*(y - 1/2) + 2000*x*y*(x - 1)*(x - 1/2)) - exp(x*y)*(2000*x*y*(y - 1)*(y - 1/2) + 2000*y*(x - 1)*(y - 1)*(y - 1/2) + 2000*y*(x - 1/2)*(y - 1)*(y - 1/2)) - x*exp(x*y)*(1000*x*(x - 1)*(x - 1/2)*(y - 1)*(y - 1/2) + 1000*x*y*(x - 1)*(x - 1/2)*(y - 1) + 1000*x*y*(x - 1)*(x - 1/2)*(y - 1/2) - 1) - y*exp(x*y)*(1000*y*(x - 1)*(x - 1/2)*(y - 1)*(y - 1/2) + 1000*x*y*(x - 1)*(y - 1)*(y - 1/2) + 1000*x*y*(x - 1/2)*(y - 1)*(y - 1/2) - 1);

mu = @(x, y) 1.0e-2;

beta = @(x, y) [1000; 1000];

sigma = @(x, y) 0.0;

f = @(x, y) - mu(x, y) * 160*tanh(40*x - 40*y)*(40*tanh(40*x - 40*y)^2 - 40) + beta(x,y)'*gradientutrue(x, y);

gD_f = @(x, y) tanh(40*(x - y))+1;

gN_f = @(x, y) mu(x, y) .* [0, 1] * gradientutrue(x, y);

%%
draw = true;
if draw
geom = defineTriangulation(K, dVertices, dBoundary, bcBoundary, ...
        bcVertices, bcValues, checkArea, checkAngle, areaValue, angleValue, true);

[A, b] = assembleSystem(K, geom, mu, beta, sigma, f, gD_f, gN_f, false);
[U, u] = computeFEsolution(geom, A, b, gD_f, false);


x = linspace(0, 1);
y = linspace(0, 1);
[X, Y] = meshgrid(x, y);
Z = utrue(X, Y);

printComparison(K, geom, U, X, Y, Z)
end

%%
areaValue = .009*[0.04, 0.01, 0.008];
errorVector = zeros(6, length(areaValue));

for exp = 1:length(areaValue)
    geom = defineTriangulation(K, dVertices, dBoundary, bcBoundary, ...
        bcVertices, bcValues, checkArea, checkAngle, areaValue(exp), angleValue, false);

    [A, b] = assembleSystem(K, geom, mu, beta, sigma, f, gD_f, gN_f);
    [U, u] = computeFEsolution(geom, A, b, gD_f);

    [err_L2, err_H1, err_Linf] = computeError(K, geom, U, utrue, gradientutrue);

    errorVector(1, exp) = sqrt(max([geom.support.TInfo(:).Area]));
    errorVector(2, exp) = err_L2;
    errorVector(3, exp) = err_H1;
    errorVector(4, exp) = err_Linf;
    errorVector(5, exp) = sqrt(min([geom.support.TInfo(:).Area]));
    errorVector(6, exp) = condest(A);
end



figure
tiledlayout(2, 3)
nexttile
semilogy(errorVector(1, :), errorVector(2, :), "b-diamond")
title("Convergenza in L^2")
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

clc
l2 = polyfit(log(errorVector(1, :)), log(errorVector(2, :)), 1);
l2 = l2(1);
h1 = polyfit(log(errorVector(1, :)), log(errorVector(3, :)), 1);
h1 = h1(1);
linf = polyfit(log(errorVector(1, :)), log(errorVector(4, :)), 1);
linf = linf(1);
cond = polyfit(log(errorVector(5, :)), log(errorVector(6, :)), 1);
cond = cond(1);
fprintf("------------------------------- \n")
fprintf("Ordine di convergenza in L2 = %d \n", round(l2))
fprintf("Ordine di convergenza in H1 = %d \n", round(h1))
fprintf("Ordine di convergenza in LInf = %d \n", round(linf))
fprintf("Ordine di crescita di cond(A) = %d \n", round(cond))
fprintf("------------------------------- \n")

