close all
clear all
clc

N = 3; 
tint = [0 10];
y0(1:N+1)=330;
y0(N+2:2*N+1)=250; 								% initial conditions

options = odeset('RelTol',1e-9,'AbsTol',[1e-9*ones(1,2*N+1)]);

[t, y] = ode45('heateq', tint, y0, options);

hold on
plot(t,y(:,2),'-r'); xlabel('Time'); ylabel('Temp');
plot(t,y(:,3),'--r'); xlabel('Time'); ylabel('Temp');
plot(t,y(:,4),'-.r'); xlabel('Time'); ylabel('Temp');
plot(t,y(:,5),'-.k'); xlabel('Time'); ylabel('Temp');
plot(t,y(:,6),'--k'); xlabel('Time'); ylabel('Temp');
plot(t,y(:,7),'-k'); xlabel('Time'); ylabel('Temp');
legend('Tube 1','Tube 2','Tube 3','Shell 3','Shell 2','Shell 1')
fig = gcf;
fig.PaperPositionMode = 'auto';
print('plot3','-dpng','-r0')


%plot(t,y(:,2),'-r'); xlabel('Time'); ylabel('Temp');
%plot(t,y(:,3),'--r'); xlabel('Time'); ylabel('Temp');
%plot(t,y(:,4),'-k'); xlabel('Time'); ylabel('Temp');
%plot(t,y(:,5),'--k'); xlabel('Time'); ylabel('Temp');
%legend('Tube 1','Tube 2','Shell 2','Shell 1')
%fig = gcf;
%fig.PaperPositionMode = 'auto';
%print('plot2','-dpng','-r0')


