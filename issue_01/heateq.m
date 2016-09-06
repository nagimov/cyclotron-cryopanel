function dy = heateq(t, y) % so y is solution function  

Tair = 298;
z = 1;
k = 450000;

% Shell Side
Cps = 4185;
T0s = 250;
As = 0.02543;
ps = 1000;
Fs = 0.1;

% Tube Side
Cpt = 1200;
T0t = 330;
At = 0.0314;
pt = 1030;
Ft = 0.2;

mt = pt * At * z;
ms = ps * As * z;

dy = zeros(2, 1);  
N  =  3;
y(1) = T0t;
y(2 * N + 2) = T0s;

for i = 2 : N + 1
    dy(i) = (pt * Ft * y(i - 1) - pt * Ft * y(i) - k * At / z * (y(i) - y(i + N + 1))) / (Cpt * mt);
end

for i = N + 2 : 2 * N + 1
    dy(i) = (ps * Fs * y(i) - ps * Fs * y(i + 1) - k * As / z * (y(i) - y(i - N - 1))) / (Cps * ms);
end

return
