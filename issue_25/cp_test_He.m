
% Cool prop wrapper comparison for Helium 

clc
close all

fluid_a = 'helium';
p = 101325 * ones(1, 100);
T = linspace(3.5, 100, 100);
h = props_htp(T, p, fluid_a, 'CP');

a = props_thp(h, p, fluid_a, 'CP')';
b = propsc_thp(h, p, fluid_a, 'CP');
delta = a-b;

subplot(2,2,1)
plot(T, delta, 'k', 'LineWidth', 2)
xlabel('Temperature of He')
ylabel('Delta T')
title('Testing THP')


c = props_uhp(h, p, fluid_a, 'CP')';
d = propsc_uhp(h, p, fluid_a, 'CP');
delta2 = c-d;

subplot(2,2,2)
plot(T, delta2, 'k', 'LineWidth', 2)
xlabel('Temperature of He')
ylabel('Delta u')
title('Testing UHP')


e = props_htp(T, p, fluid_a, 'CP');
f = propsc_htp(T, p, fluid_a, 'CP');
delta2 = e-f;

subplot(2,2,3)
plot(T, delta2, 'k', 'LineWidth', 2)
xlabel('Temperature of He')
ylabel('Delta h')
title('Testing HTP')
