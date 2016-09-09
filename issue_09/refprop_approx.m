clc
clear all

%Pressure and quality
q_a = .75;
p_b = 101; 

% Generate a range of temperatures for nitrogen and helium
T_a = linspace(75,80,100);
T_b = linspace(4.3,6,100);

% Compute h and u using refprop  
for i=1:100
h_a(i) = refpropm('U','T',T_a(i),'Q',q_a,'nitrogen'); 
h_b(i) = refpropm('U','T',T_b(i),'P',p_b,'helium'); 
u_a(i) = refpropm('U','T',T_a(i),'Q',q_a,'nitrogen');
u_b(i) = refpropm('U','T',T_b(i),'P',p_b,'helium');
end

% Fit a line using 'polyfit' 
% In each case compute r2 values to see how good the fit was 
p = polyfit(T_a,h_a,1);
yfit = polyval(p,T_a);
yresid = h_a - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(h_a)-1) * var(h_a);
rsq(1) = 1 - SSresid/SStotal;

p = polyfit(T_a,u_a,1);
yfit = polyval(p,T_a);
yresid = u_a - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(u_a)-1) * var(u_a);
rsq(2) = 1 - SSresid/SStotal;

p = polyfit(T_b,h_b,1);
yfit = polyval(p,T_b);
yresid = h_b - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(h_b)-1) * var(h_b);
rsq(3) = 1 - SSresid/SStotal;

p = polyfit(T_b,u_b,1);
yfit = polyval(p,T_b);
yresid = u_b - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(u_b)-1) * var(u_b);
rsq(4) = 1 - SSresid/SStotal;

%Plot 
figure
subplot(2,2,1)
plot(T_a,h_a,'-b','LineWidth',2); xlabel('Temp'); ylabel('h_a');
subplot(2,2,2)
plot(T_b,h_b,'-b','LineWidth',2); xlabel('Temp'); ylabel('h_b');
subplot(2,2,3)
plot(T_a,u_a,'-b','LineWidth',2); xlabel('Temp'); ylabel('u_a');
subplot(2,2,4)
plot(T_b,u_b,'-b','LineWidth',2); xlabel('Temp'); ylabel('u_a');
fig = gcf;
fig.PaperPositionMode = 'auto';
print('refprop_approx','-dpng','-r0')
close all 

