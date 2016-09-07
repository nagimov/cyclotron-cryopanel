close all
clear all
clc
global mid last 

%How many stages?
N = 3; 
% Compute key indicies 
mid = N + 2;               %where He columns start
last = 2 * (N +  1);       %last column

% What is our time interval?
t_interval = [0 10];
% Initial conditions 
T0(1 : mid - 1) = 80;
T0(mid : last) = 4; 							

% Some technical options for ode45 - to ensure accuracy 
options = odeset('RelTol',1e-3,'AbsTol',[1e-3*ones(1,2*( N + 1 ))]);
% Solve 
% The solver produces matrix T with solution data. 
% Each row contains T at each stage of the heat exchanger, for a fixed t.
% For N=3, columns 1-4 contain stream A, columns 5-8 contain stream B. 
[t, T] = ode45('heateq', t_interval, T0, options);

% ******************
% Plot T vs time 
hold on
plot(t,T(:,2),'-r','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(t,T(:,3),'--r','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(t,T(:,4),'-.r','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(t,T(:,6),'-.b','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(t,T(:,7),'--b','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(t,T(:,8),'-b','LineWidth',2); xlabel('Time'); ylabel('Temp');
legend('N2 Stage 1','N2 Stage 2','N2 Stage 3','He Stage 3','He Stage 2','He Stage 1','location','eastoutside')
fig = gcf;
fig.PaperPositionMode = 'auto';
print('plot_T_vs_time','-dpng','-r0')

% ******************
% Plot T vs length of heat exchanger stages  
figure 
hold on
a = round(length(t)/2); %find the middle of time vector
plot(1:4,T(a,1:4),'-r','LineWidth',2); xlabel('Stage'); ylabel('Temp');
plot(1:4,fliplr(T(a,5:8)),'-b','LineWidth',2); xlabel('Stage'); ylabel('Temp');
%we need to flip the He vector since we are counter-current
legend('Nitrogen','Helium','location','eastoutside')
fig = gcf;
fig.PaperPositionMode = 'auto';
print('plot_T_vs_length','-dpng','-r0')
