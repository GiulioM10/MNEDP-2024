function [t, Dt] = discretizeTime(t0, T, approxDt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
N = ceil((T-t0)/approxDt);
t = linspace(t0, T, N+1);
Dt = (T-t0)/N;
end