clc
clear all
global t mid last 

t = 5; %How many time steps do we take? (time step be specified, currently 0.1s)
N = 4; %How many slices?

%Calculate indicies
mid = N + 2;               %where He columns start
last = 2 * (N +  1);       %last column

%Take an intial guess 
T0a = 80 * ones ( t + 1 , mid - 1 );
T0b = 4.5 * ones ( t + 1 , mid - 1 ); 
T0 = [T0a T0b];

%Constraints on temperature - this ensures that temperatures iterated are
%always within refprop limits!!! 
Tmax_a = 90 * ones ( t + 1 , mid - 1 );
Tmax_b = 5 * ones ( t + 1 , mid - 1 );
Tmax = [Tmax_a Tmax_b];
Tmin_a = 70 * ones ( t + 1 , mid - 1 );
Tmin_b = 3 * ones ( t + 1 , mid - 1 );
Tmin = [Tmin_a Tmin_b];

%Constraints on the intial data - this ensures that values are always fix
%at 80 and 4.5!!! 
Tmax(1, 1 : mid - 1) = 80;
Tmax(1, mid : last) = 4.5; 
Tmax(1 : t + 1, 1) = 80;
Tmax(1 : t + 1, mid) = 4.5;
Tmin(1, 1 : mid - 1) = 80;
Tmin(1, mid : last) = 4.5; 
Tmin(1 : t + 1, 1) = 80;
Tmin(1 : t + 1, mid) = 4.5;

%Solve 
options = optimset('TolX', 1e-6, 'TolFun', 1e-6, ...
    'MaxFunEvals', 1e7, 'MaxIter', 1e7);
[T,fval,exitflag,output] = fmincon ( 'heateq' , T0  , [], [], [], [], Tmin, Tmax, [] , options); 
% If we solve taking into account constraints, we must use fmincon. There
% are a lot of paramaters for fmincon that we do not use, so we just insert
% place holders [] into those spots

% Old code for fsolve
% [T,fval,exitflag,output] = fsolve ( 'heateq' , T0 , options);


%Plot
hold on
plot(1:t+1,T(:,2),'-r','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(1:t+1,T(:,3),'--r','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(1:t+1,T(:,4),'-.r','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(1:t+1,T(:,6),'-.b','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(1:t+1,T(:,7),'--b','LineWidth',2); xlabel('Time'); ylabel('Temp');
plot(1:t+1,T(:,8),'-b','LineWidth',2); xlabel('Time'); ylabel('Temp');
legend('N2 Slice 1','N2 Slice 2','N2 Slice 3','He Slice 3','He Slice 2','He Slice 1','location','eastoutside')
fig = gcf;
fig.PaperPositionMode = 'auto';
print('plot_T_vs_time','-dpng','-r0')

figure 
hold on
plot(1:mid-1,T(6,1:mid-1),'-r','LineWidth',2); xlabel('Slice'); ylabel('Temp');
plot(1:mid-1,fliplr(T(6,mid:last)),'-b','LineWidth',2); xlabel('Slice'); ylabel('Temp');
%we need to flip the He vector since we are counter-current
legend('Nitrogen','Helium','location','eastoutside')
fig = gcf;
fig.PaperPositionMode = 'auto';
print('plot_T_vs_length','-dpng','-r0')
