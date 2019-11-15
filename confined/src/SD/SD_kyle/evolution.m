function [tout, rout, qout] = evolution(t)
%EVOLUTION Integrates equation of motion for charged particle in a channel.
%   Calculates trajectory of a single, conductive sphere suspended in a
%   viscous fluid between two parallel, oppositely charged plates.

global H      % channel height
global L      % channel width
global D      % channel depth
global Pe     % Peclet number
global drift  % drift velocity function

H = 3;
L = 10;
D = L;
Pe = 1;

% 
qel = [-1, 1];

delta = 1.0e-02;

% Initialize particle position and charge.
[r, q] = setup(L, H, D);

% Initialize matrices for storing particle trajectory.
tout = [ 0 ];
rout = [ r' ];
qout = [ q ];

% For a single particle, drift velocity does not change in time as it depends
% only on particle position and channel geometry (heigth), thus its
% functional from can be pre-computed outside the integration loop.
dh = 0.005;
[h, v] = vdrift(H, dh);
drift = griddedInterpolant(h, v);

% Integerate equation of motion over a given time interval [0; tfin].
dt = 0.01;
dtau = 0;
tfin = t;
t = 0;
while t < tfin;
    % Calculate mobility matrix.
    M = mobility(r);
    
    % Find out a matrix $A$ such as $A A^{T} = M$.
    A = chol(M);
    
    % Calculate new position of the particle.
    v = Pe * M * force(r, q);
    vd = drift(r(3)) * [0; 0; 1];
%    rrnd = mvnrnd(zeros(3,1), 2 * M)';
    rrnd = sqrt(2) * A * randn(3, 1);
        
    dr = v * dt + vd * dt + rrnd * sqrt(dt);
    rnew = r + dr;
    
    % Check if any events (wall crossing, periodic boundary crossing) occured
    % during the time step and take appropriate actions.
    value = events(rnew);
    ie = find(value);
    for i = 1:length(ie)
        % Crossing any of the electrodes.
        if ie(i) < 3
            % Find out smaller time step which leads to a contact of the
            % sphere with an electrode.
            dhnew = norm(rnew(3) - r(3));
            if ie == 1
                dhcol = norm((1 + delta) - r(3));
            else
                dhcol = norm(H - (1 + delta) - r(3));
            end
            x = roots([norm(v + vd), norm(rrnd), -dhcol / dhnew * norm(dr)]);

            % Find the smallest, positive root.
            dtau = min(x(x > 0))^2;

            % Recalculate new position using smaller time step.
            rnew = r + v * dtau + vd * dtau + rrnd * sqrt(dtau);
           
            % Update charge.
            q = qel(ie(i));
            
        % Crossing periodic boundries.
        else
            rnew = rnew - [round(r(1) / L) * L; round(r(2) / D) * D; 0];
        end
    end
    
    % Advance clock.
    if dtau == 0
        t = t + dt;
    else
        t = t + dtau;
        dtau = 0;
    end
    r = rnew;
    
    % Accumulate trajectory.
    tout = [tout; t];
    rout = [rout; r'];
    qout = [qout; q];
end

end


function [value] = events(r)
%EVENTS Defines possible events which may occur in the system.
%   Locates times 

global D H L;

delta = 1.0e-02;

% Initialize event arrays:
%
% (a) collision with lower and upper electrode, respectively;
value = [r(3) - 1 < delta; (H - 1) - r(3) < delta; ...
         r(1) < delta; r(1) > L - delta; ...
         r(2) < delta; r(2) > D - delta];
end