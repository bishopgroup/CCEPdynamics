function [f] = force(r, q)
%FORCE Calculates an external, deterministic force acting on a particle.
%   Returns force acting on a particle due to the presence of the field
%   between the electrodes.

f = q * [0; 0; -1];
end

