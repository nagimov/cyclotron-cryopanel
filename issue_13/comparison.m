clear 
clc

HX_slices = 10;

% Compute T without time-step adjustment 
[T_a_sol4, T_b_sol4, T_w_sol4] = RN_04;
close all 

% Compute T with time-step adjustment 
[T_a_sol5, T_b_sol5, T_w_sol5] = RN_05;
close all

% Compute differences in temperatures 
deltaT_a = T_a_sol4 - T_a_sol5;
deltaT_b = T_b_sol4 - T_b_sol5;
deltaT_w = T_w_sol4 - T_w_sol5; 

% PLOT 
figure
hold on
plot(1:HX_slices,deltaT_a,'r')
plot(1:HX_slices,deltaT_b,'b')
plot(1:HX_slices,deltaT_w,'g')
%axis([1 HX_slices -15 15])
xlabel('Slices')
ylabel('Temperature Delta')
fig = gcf;
fig.PaperPositionMode = 'auto';
print(['plot_deltaT' num2str(HX_slices)],'-dpng','-r0')
