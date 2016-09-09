function F = heateq(T)  %Here F compute difference between our equation and zero i.e. F(T)~0. 
global t mid last 

%Time Step
delta_t = .001;

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
q_a = .75;
p_a = 101; 
k_a = 1e6;
k_b = 1e6; 

%Get U for this data 
U_a(1,1) = refpropm('U','T',T(1,1),'Q',q_a,'nitrogen');
U_b(1,1) = refpropm('U','T',T(1,mid),'P',p_a,'helium');

%Equations Nitrogen 
for j = 2 : mid - 1 %slices
for i = 2 : t + 1 %time step 

delta_T = T(i, last - j + 1 ) - T(i,j);
U_a(i,j) = refpropm('U','T',T(i,j),'Q',q_a,'nitrogen');
H_a_in(i,j) = refpropm('H','T',T(i,j-1),'Q',q_a,'nitrogen');
H_a_out(i,j) = refpropm('H','T',T(i,j),'Q',q_a,'nitrogen');
Q_cond_a = k_a * As / HX_L * delta_T ;

F(i,j) = U_a(i-1,j) - U_a(i,j) + delta_t / M * ( m * H_a_in(i,j) - m * H_a_out(i,j) + Q_cond_a ); 

end 
end 

%Equations Helium  
for j = mid + 1 : last %slices
for i = 2 : t + 1  %time step 

delta_T = T(i, last - j + 1 ) - T(i,j); 
U_b(i,j) = refpropm('U','T',T(i,j),'P',p_a,'helium');
H_b_in(i,j) = refpropm('H','T',T(i,j-1),'P',p_a,'helium');
H_b_out(i,j) = refpropm('H','T',T(i,j),'P',p_a,'helium');
Q_cond_b = k_b * As / HX_L * delta_T;

F(i,j) = U_b(i-1,j) - U_b(i,j) + delta_t / M * ( m * H_b_in(i,j) - m * H_b_out(i,j) + Q_cond_b ); 

end 
end

%Ensure boundary conditions
F(1, 1 : mid - 1) = 80 - T(1, 1 : mid - 1);
F(1, mid : last) = 4.5 - T(1, mid : last); 
F(t + 1, 1) = 80 - T(t + 1 ,1);
F(t + 1, mid) = 4.5 - T(i, mid);
