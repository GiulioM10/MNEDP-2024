clear all
close all
clc

syms x y
N1(x, y) = x;
N2(x, y) = y;
N3(x, y) = 1 - x -y;

Phi1(x, y) = 2 * N1(x, y) * (N1(x, y) - 0.5);
Phi2(x, y) = 2 * N2(x, y) * (N2(x, y) - 0.5);
Phi3(x, y) = 2 * N3(x, y) * (N3(x, y) - 0.5);
Phi4(x, y) = 4 * N1(x, y) * N3(x, y);
Phi5(x, y) = 4 * N1(x, y) * N2(x, y);
Phi6(x, y) = 4 * N2(x, y) * N3(x, y);