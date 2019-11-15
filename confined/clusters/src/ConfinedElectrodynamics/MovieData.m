load('RESULTS.mat')

t = linspace(t_Save(1),t_Save(end),length(t_Save))';
x123 = interp1(t_Save,x123_Save,t,'spline');
q = interp1(t_Save,sum(q_Save,2),t,'nearest');
q = [q,q,q];

output = [t,q,x123];
save('TrimerData.txt','output','-ascii','-double','-tabs');

% for i = 1:size(x123,1)
%     draw(reshape(x123(i,:),3,part.Np), q(i,:), param);
%     pause(0.01);
% end

x = mean(x123(:,1),2);
plot(x,t,'-o') 
