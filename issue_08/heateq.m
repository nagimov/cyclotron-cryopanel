function F = heateq(T)  %Here F compute difference between our equation and zero i.e. F(T)~0. 
global t mid last 

%Time Step
delta_t = .1;

%Data
m = 100; 
M = 100; 

% Area Information
diam = 0.0254;
rad = diam/2;
HX_L = 1;
Ac = pi * rad^2;
As = 2 * pi * rad * HX_L;

% Heat Transfer Coefficients 
k_a = 1e6;
k_b = 1e6; 

%Get U for this data 
U_a(1,1) = refpropm('U','T',T(1,1),'Q',.75,'nitrogen');
U_b(1,1) = refpropm('U','T',T(1,mid),'P',101,'helium');

%Equations Nitrogen 
for j = 2 : mid - 1 %stages
for i = 2 : t + 1 %time step 

delta_T = T(i, last - j + 1 ) - T(i,j);
U_a(i,j) = refpropm('U','T',T(i,j),'Q',.75,'nitrogen');
H_a_in(i,j) = refpropm('H','T',T(i,j-1),'Q',.75,'nitrogen');
H_a_out(i,j) = refpropm('H','T',T(i,j),'Q',.75,'nitrogen');
Q_cond_a = k_a * As / HX_L * delta_T ;

F(i,j) = U_a(i-1,j) - U_a(i,j) + delta_t * ( m * H_a_in(i,j) - m * H_a_out(i,j) + Q_cond_a ) / M; 

end 
end 

%Equations Helium  
for j = mid + 1 : last %stages
for i = 2 : t + 1  %time step 

delta_T = T(i, last - j + 1 ) - T(i,j); 
U_b(i,j) = refpropm('U','T',T(i,j),'P',101,'helium');
H_b_in(i,j) = refpropm('H','T',T(i,j-1),'P',101,'helium');
H_b_out(i,j) = refpropm('H','T',T(i,j),'P',101,'helium');
Q_cond_b = k_b * As / HX_L * delta_T;

F(i,j) = U_b(i-1,j) - U_b(i,j) + delta_t * ( m * H_b_in(i,j) - m * H_b_out(i,j) + Q_cond_b ) / M; 

end 
end

%Ensure boundary conditions
for j = 1 : mid - 1
F(1, j) = 1e6*(80 - T(1,j));
end

for j = mid : last 
F(1, j) = 1e6*(4 - T(1,j)); 
end

for i = 1 : t + 1
F(i, 1) = 1e6*(80 - T(i,1));
end

for i = 1 : t + 1
F(i, mid) = 1e6*(4 - T(i, mid));
end




