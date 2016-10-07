function [x,fval,exitflag] =  heateq(j) %Time iteration 

global mid last Hsol delta_t p 
options = optimset('TolX', 1e-7, 'TolFun', 1e-7, ...
    'MaxFunEvals', 1e7, 'MaxIter', 1e7, 'Display', 'iter');
x0 = Hsol(j-1,:); 

[x,fval,exitflag] = fsolve ( @eqgen , x0 , options);

%[],[],[],[],[],[],[],
function F = eqgen(H)  %Here F compute difference between our equation and zero i.e. F(H)~0. 

%Data
m = 1; 
M = 1; 

% Area Information
diam = 0.0254;
radius = diam/2;
HX_L = 1;
Ac = pi * radius^2;
As = 2 * pi * radius * HX_L;

% Thermodynamic Parameters  
p_a = 100; 
p_b = 100; 
k_a = 1e4;
k_b = 1e4; 

%Equations Nitrogen 
for i = 2 : mid - 1 %slices
   
U_a(i) = polyval(p(:,2),H(i));
U_a_ini(i) = polyval(p(:,2),Hsol(j-1,i)); 
T(i) = polyval(p(:,1),H(i)); %T for stream A
T(last - i + 1) = polyval(p(:,3),H(last - i + 1)); %T for stream B
delta_T = T(last - i + 1) - T(i);
Q_cond_a = k_a * As / HX_L * delta_T;

F(i) = U_a(i) - U_a_ini(i) + delta_t / M * ( m * H(i-1) - m * H(i) + Q_cond_a ); 
end 
 

%Equations Helium  
for i = mid + 1 : last %slices
    
U_b(i) = polyval(p(:,4),H(i));
U_b_ini(i) = polyval(p(:,4),Hsol(j-1,i)); 
T(i) = polyval(p(:,3),H(i)); %T for stream B
T(last - i + 1) = polyval(p(:,1),H(last - i + 1)); %T for stream A
delta_T = T(last - i + 1) - T(i); 
Q_cond_b = k_b * As / HX_L * delta_T;

F(i) = U_b(i) - U_b_ini(i) + delta_t / M * ( m * H(i-1) - m * H(i) + Q_cond_b ); 
end

%Boundary constraints 
A = refpropm('H','T',80,'P',p_a,'nitrogen'); %Convert H into T on the boundary 
B = refpropm('H','T',4.5,'P',p_b,'helium'); %Convert H into T on the boundary 
F(1) = A - H(1);
F(mid) = B - H(mid);

eval = sum(F);
end
end 
