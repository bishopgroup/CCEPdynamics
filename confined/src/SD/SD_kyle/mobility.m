function [ M ] = mobility(r)
%MOBILITY Calculates instantenous mobility matrix.
%   Detailed explanation goes here

global H

% Calculate fractional height, $z / H$, of the particle.
Xi = r(3) / H;

% Find out the grand resitstance matrix.
R = RSD(Xi, H);

RFU = R(1:3, 1:3);
M = inv(RFU);
end