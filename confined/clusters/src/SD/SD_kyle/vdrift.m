function [h, v] = vdrift(H, dh)
%DRIFT Computes the drift velocity of a signle particle in a channel
%
%   drift(H) computes the drift velocity of a single, Brownian particle in
%   a channel of width $H / a$ where $a$ is particle's radius.

hmin = 1;
hmax = H - 1;
h = hmin:dh:hmax;

% Discard first element as it corresponds to sphere being in contact with
% the wall.
h(1) = [];
h(end) = [];

% Express sphere positions above the wall as fractions of channel height.
xi = h / H;

n = length(xi);
rows = zeros(n, 3);
for i = 1:n
    R = RSD(xi(i), H);
    RFU = R(1:3, 1:3);
    MFU = inv(RFU);
    
    rows(i, :) = MFU(3, :);
end

h(end) = [];
h(1) = [];
v = cdiff(rows(:,3), dh);
end