clc
clear all
global mid last Tsol delta_t 

delta_t = .01; %Time Step
t = 20; %How many time steps do we take? 
N = 10; %How many slices?

%Calculate indicies
mid = N + 2;               %where He columns start
last = 2 * (N +  1);       %last column

%Take first intial guess 
T0a = 80 * ones ( 1 , mid - 1 );
T0b = 4.5 * ones ( 1 , mid - 1 ); 
Tsol(1, :) = [T0a T0b];

%Solve 
for j = 2 : t + 1
[Tsol(j, :),fval(j,:),exitflag(j,:)] = heateq(j); 
end 

hold on
plot(1:11,Tsol(:,1:mid-1))
plot(1:11,Tsol(:,mid:last))
