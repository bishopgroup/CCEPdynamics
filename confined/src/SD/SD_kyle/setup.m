function [r, q] = setup(L, H, D)
%SETUP Initializes particle's position and charge.
%   Puts particle in the center of the channel and assigns a positive
%   charge to it.

r = [L/2; D/2; H/2];
q = -1;
end

