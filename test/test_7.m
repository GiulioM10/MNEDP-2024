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

K = 1;

%%
utrue_t = @(x, y, t) sin(pi * t/2) + x + y;

dt_utrue_t = @(x, y, t) (pi*cos((pi*t)/2))/2;

gradientutrue_t = @(x, y, t) [1;
                              1];

laplacianutrue = @(x, y, t) 2;

rho = @(x, y) cosh(x - y);

mu = @(x, y) x + y + 1;

beta = @(x, y) [-1; x*y];

sigma = @(x, y) 5;

f = @(x, y, t) rho(x, y)*dt_utrue_t(x, y, t) -laplacianutrue(x, y, t) + beta(x, y)'*gradientutrue_t(x,y,t) + sigma(x, y)*utrue_t(x, y, t);

gD_f = @(x, y, t) sin(pi * t/2) + x + y;

Dt_gD_f = @(x, y, t) (pi*cos((pi*t)/2))/2;

gN_f = @(x, y, t) mu(x, y) .* [0, 1] * gradientutrue_t(x, y, t);

u0_f = @(x, y) x + y;

%%
approxDt =  [0.1, 0.05, 0.01, 0.005, 0.001];
times = zeros(2, length(approxDt));
geom = defineTriangulation(K, dVertices, dBoundary, bcBoundary, ...
    bcVertices, bcValues, checkArea, checkAngle, .04*0.04, angleValue, false);
for type = 1:2
for exp = 1:length(approxDt)

    [t, dt] = discretizeTime(t0, T, approxDt(exp));
    tic
    [~, ~, ~] = solveEvoProblem(K, geom, t, dt, rho, mu, beta, sigma, f, gD_f, Dt_gD_f, gN_f, u0_f, false, false, ~(type==1));
    time = toc;

    times(type, exp) = time;
end
end


cumul_times = cumsum(times, 2);

X = categorical({'0.1','0.05','0.01','0.005', '0.001'});
X = reordercats(X,{'0.1','0.05','0.01','0.005', '0.001'});

figure
bar(X, log(1+times'))
legend("Direct solver", "Iterative solver", 'Location', 'northwest')
xlabel("\Delta t")
ylabel("log(1 +time)")

% figure
% plot(1:length(approxDt), cumul_times(1,:), "-x", 1:length(approxDt), cumul_times(2,:), "-+")
% grid on
% legend("Direct solver", "Iterative solver", 'Location', 'northwest')
% ylabel("time [s]")

