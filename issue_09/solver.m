tic
clc
clear all
close all 
global mid last Hsol delta_t p

delta_t = .01; %Time Step
t = 10; %How many time steps do we take? 
N = 5; %How many slices?

%Calculate indicies
mid = N + 2;               %where He columns start
last = 2 * (N +  1);       %last column

%Take first intial guess 
p_a = 100; %Nitrogen Pressure 
p_b = 100; %Helium Pressure 
A = refpropm('H','T',80,'P',p_a,'nitrogen'); %To get an initial guess for H, since we know T 
B = refpropm('H','T',4.5,'P',p_b,'helium'); %To get an initial guess for H, since we know T 
H0a = A * ones ( 1 , mid - 1 );
H0b = B * ones ( 1 , mid - 1 ); 
Hsol(1, :) = [H0a H0b];

%Linear refprop
p = refprop_approx(p_a,p_b); 

%Solve 
for j = 2 : t + 1
[Hsol(j, :),fval(j,:),exitflag(j,:)] = heateq(j); 
disp(['Time step completed ' num2str(j)])
end 

%Plot 1: Temp vs. length  
hold on
plot(1:N+1,Hsol(:,1:mid-1))
plot(1:N+1,Hsol(:,fliplr(mid:last)))
xlabel('Slices')
ylabel('Enthalpy')
fig = gcf;
fig.PaperPositionMode = 'auto';
print('plot_T_vs_length','-dpng','-r0')

%Plot 2: Temp vs. time 
figure
hold on 
plot(1:t,Hsol(1:t,1:mid-1))
plot(1:t,Hsol(1:t,mid:last))
xlabel('Time')
ylabel('Enthalpy')
fig = gcf;
fig.PaperPositionMode = 'auto';
print('plot_T_vs_time','-dpng','-r0')
toc 
