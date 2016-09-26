function [x,fval,exitflag] =  heateq(j) 

global mid last Tsol delta_t 
options = optimset('TolX', 1e-7, 'TolFun', 1e-7, ...
    'MaxFunEvals', 1e7, 'MaxIter', 1e7, 'Display', 'off');
x0 = Tsol(j-1,:); 

[x,fval,exitflag] = fsolve ( @eqgen , x0  , options);


function F = eqgen(T)  %Here F compute difference between our equation and zero i.e. F(T)~0. 

%Data
m = 1000; 
M = 1; 

% Area Information
diam = 0.0254;
radius = diam/2;
HX_L = 1;
Ac = pi * radius^2;
As = 2 * pi * radius * HX_L;

% Thermodynamic Parameters  
q_a = .75;
p_a = 10; 
k_a = 1e5;
k_b = 1e5; 

%Equations Nitrogen 
for i = 2 : mid - 1 %slices

delta_T = T(last - i + 1 ) - T(i);
U_a(i) = refpropm('U','T',T(i),'P',p_a,'nitrogen');
U_a_ini(i) = refpropm('U','T',Tsol(j-1,i),'P',p_a,'nitrogen');
H_a_in(i) = refpropm('H','T',T(i-1),'P',p_a,'nitrogen');
H_a_out(i) = refpropm('H','T',T(i),'P',p_a,'nitrogen');
Q_cond_a = k_a * As / HX_L * delta_T ;

F(i) = U_a(i) - U_a_ini(i) + delta_t / M * ( m * H_a_in(i) - m * H_a_out(i) + Q_cond_a ); 

end 
 

%Equations Helium  
for i = mid + 1 : last %slices
    
delta_T = T(last - i + 1) - T(i); 
U_b(i) = refpropm('U','T',T(i),'P',p_a,'helium');
U_b_ini(i) = refpropm('U','T',Tsol(j-1,i),'P',p_a,'helium');
H_b_in(i) = refpropm('H','T',T(i-1),'P',p_a,'helium');
H_b_out(i) = refpropm('H','T',T(i),'P',p_a,'helium');
Q_cond_b = k_b * As / HX_L * delta_T;

F(i) = U_b(i) - U_b_ini(i) + delta_t / M * ( m * H_b_in(i) - m * H_b_out(i) + Q_cond_b ); 

end

F(1) = 80 - T(1);
F(mid) = 4.5 - T(mid);

end
end 
