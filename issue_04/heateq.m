function dT = heateq(t, T) 
global N 

% A - Nitrogen Data
T_a_in = 80;
rho_a = 1000;
F_a = 2e-4;

% B - Helium Data
T_b_in = 4;
rho_b = 1030;
F_b = 2e-4;

% Area Information
diam = 0.0254;
Ac = pi*(diam/2)^2;
HX_L = 1;

% Heat Transfer Coefficients 
HX_UA_a = 2100;
HX_UA_b = 2100; 

% Compute Masses
m_a = rho_a*Ac*HX_L;
m_b = rho_b*Ac*HX_L;

% Prepare to solve - preallocate dummy vector
dT = zeros(2, 1);  
% Prepare to solve - enter initial data into the DE 
T(1) = T_a_in;
T(N + 2) = T_b_in;

% A - Nitrogen DE
for i = 2 : N + 1
    Cp_a = refpropm('C','T',T(i),'Q',0,'nitrogen');
    Q_in_a = rho_a * F_a * Cp_a * T(i - 1); 
    Q_out_a = rho_a * F_a * Cp_a * T(i);
    dT(i) = (Q_in_a - Q_out_a - HX_UA_a) / (Cp_a * m_a);
end

% B - Helium DE
for i = N + 3 : 2 * N + 2
    Cp_b = refpropm('C','T',T(i),'Q',0,'helium'); 
    Q_in_b = rho_b * F_b * Cp_b * T(i - 1);
    Q_out_b = rho_b * F_b * Cp_b *  T(i);
    dT(i) = (Q_in_b - Q_out_b + HX_UA_b) / (Cp_b * m_b);
end

return
