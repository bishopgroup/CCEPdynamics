function [frac,flag,adjmat] = Overlap(Parameter,partNEXT,Particle)
%{
%}

flag = false;
frac = 1;
adjmat = Particle.adjmat;

for i = 1:partNEXT.nParticle
    %% Loop through particle pairs
%     for j = i+1:partNEXT.nParticle
%         % Interparticle separation (r = xj-xi)
%         dx1 = partNEXT.position(1,j)-partNEXT.position(1,i);
%         dx2 = partNEXT.position(2,j)-partNEXT.position(2,i);
%         dx3 = partNEXT.position(3,j)-partNEXT.position(3,i);
%         dx1 = dx1 - round(dx1 / Parameter.L1) * Parameter.L1; % periodic bc
%         dx2 = dx2 - round(dx2 / Parameter.L2) * Parameter.L2; % periodic bc
%         e = [dx1,dx2,dx3];
%         rNEXT = norm(e);
%         
%         if rNEXT < 2+Parameter.delta % there is overlap
%             flag = true;
%             adjmat(i,j) = 1; 
%             adjmat(j,i) = 1;
%         
%             % Compute Interparticle separation at previous time step (r = xj-xi)
%             dx1 = Particle.position(1,j)-Particle.position(1,i);
%             dx2 = Particle.position(2,j)-Particle.position(2,i);
%             dx3 = Particle.position(3,j)-Particle.position(3,i);
%             dx1 = dx1 - round(dx1 / Parameter.L1) * Parameter.L1; % periodic bc
%             dx2 = dx2 - round(dx2 / Parameter.L2) * Parameter.L2; % periodic bc
%             e = [dx1,dx2,dx3];
%             r = norm(e);
%             
%             tmp = (r - (2+Parameter.delta)) / (r - rNEXT);
%             if tmp < frac
%                 frac = tmp;
%             end
%         end
%     end
    
    %% Wall Interactions with Lower Wall (image approach)
    dx3 = -2*partNEXT.position(3,i);
    e = [0,0,dx3];
    rNEXT = abs(dx3);
    
    if round(rNEXT,2) < 2+Parameter.delta % there is overlap
        flag = true;
        adjmat(i,i) = 2; 
        
        % Compute Interparticle separation at previous time step (r = xj-xi)
        dx3 = -2*Particle.position(3,i);
        e = [0,0,dx3];
        r = abs(dx3);
        
        tmp = (r - (2+Parameter.delta)) / (r - rNEXT);
        if tmp < frac
            frac = tmp;
        end
    end

    
    %% Wall Interactions with Upper Wall (image approach)
    dx3 = 2*(Parameter.L3-partNEXT.position(3,i));
    e = [0,0,dx3];
    rNEXT = abs(dx3);

    if round(rNEXT,2) < 2+Parameter.delta % there is overlap
        flag = true;
        adjmat(i,i) = 3;
           
        % Compute Interparticle separation at previous time step (r = xj-xi)
        dx3 = 2*(Parameter.L3-Particle.position(3,i));
        e = [0,0,dx3];
        r = abs(dx3);
        
        tmp = (r - (2+Parameter.delta)) / (r - rNEXT);
        if tmp < frac
            frac = tmp;
        end
    end
end
