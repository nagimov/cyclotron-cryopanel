function F = heateq(T)  %Here F compute difference between our equation and zero i.e. F(T)~0. 
global t N

%Time Step
delta_t = 0.1; 

%Data
m = 1000; 
HX_UA_a = 6e6; %For now assume fixed flux
HX_UA_b = 6e6; %For now assume fixed flux

%Set initial data
T(1, 1 : N + 1 ) = 80; 
T(1, N + 2 : 2 * ( N + 1 )) = 4.5; 
T(:, 1) = 80;
T(:, N + 2) = 4.5; 
%Get U for this data 
U_a(1) = refpropm('U','T',T(1,1),'P',101,'nitrogen');
U_b(1) = refpropm('U','T',T(1,N+2),'P',101,'helium');

%Equations Nitrogen 
for j = 2 : N + 1 %stages
for i = 1 : t + 1 %time step 
    
U_a(i+1) = refpropm('U','T',T(i,j),'P',101,'nitrogen');
H_a_in = refpropm('H','T',T(i,j-1),'P',101,'nitrogen');
H_a_out = refpropm('H','T',T(i,j),'P',101,'nitrogen');

F(i+1,j) = U_a(i+1) - U_a(i) + delta_t * ( m * H_a_in - m * H_a_out - HX_UA_a ); 

end 
end 

%Equations Helium  
for j = N + 3 : 2*(N + 1) %stages
for i = 1 : t + 1  %time step 
    
U_b(i+1) = refpropm('U','T',T(i,j),'P',101,'helium');
H_b_in = refpropm('H','T',T(i,j-1),'P',101,'helium');
H_b_out = refpropm('H','T',T(i,j),'P',101,'helium');

F(i+1,j) = U_b(i+1) - U_b(i) + delta_t * ( m * H_b_in - m * H_b_out + HX_UA_b ); 

end 
end

