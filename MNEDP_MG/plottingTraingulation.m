function [T2] = plottingTraingulation(T)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ntriang = size(T, 1);
T2 = zeros(4*ntriang, 3);
count = 0;

for i = 1:ntriang
    count = count + 1;
    T2(count, :) = [T(i, 1), T(i, 5), T(i, 4)];
    count = count + 1;
    T2(count, :) = [T(i, 4), T(i, 5), T(i, 6)];
    count = count + 1;
    T2(count, :) = [T(i, 4), T(i, 6), T(i, 3)];
    count = count + 1;
    T2(count, :) = [T(i, 5), T(i, 2), T(i, 6)];
end
end