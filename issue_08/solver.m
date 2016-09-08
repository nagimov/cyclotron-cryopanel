clc
clear all
global t N 

t = 5; %How many time steps do we take? (time step be specified, currently 0.1s)
N = 3; %How many stages?

%Take an intial guess 
T0a = 80 * ones ( t + 1 , N + 1 );
T0b = 4.5 * ones ( t + 1 , N + 1 ); 
T0 = [T0a T0b ];

%Solve 
T = fsolve ( 'heateq' , T0 ); 

%Plot
hold on
plot(1:t+1,T(:,2),'-r','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(1:t+1,T(:,3),'--r','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(1:t+1,T(:,4),'-.r','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(1:t+1,T(:,6),'-.b','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(1:t+1,T(:,7),'--b','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(1:t+1,T(:,8),'-b','LineWidth',2); xlabel('Time'); ylabel('Temp');
legend('N2 Stage 1','N2 Stage 2','N2 Stage 3','He Stage 3','He Stage 2','He Stage 1','location','eastoutside')
fig = gcf;
fig.PaperPositionMode = 'auto';
print('plot_T_vs_time','-dpng','-r0')

figure 
hold on
plot(1:N+1,T(6,1:N+1),'-r','LineWidth',2); xlabel('Stage'); ylabel('Temp');
plot(1:N+1,fliplr(T(6,N+2:2*(N+1))),'-b','LineWidth',2); xlabel('Stage'); ylabel('Temp');
%we need to flip the He vector since we are counter-current
legend('Nitrogen','Helium','location','eastoutside')
fig = gcf;
fig.PaperPositionMode = 'auto';
print('plot_T_vs_length','-dpng','-r0')
