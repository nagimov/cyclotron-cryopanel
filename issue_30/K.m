function eval = K(T)  % W/(m*K)

% http://cryogenics.nist.gov/MPropsMAY/304Stainless/304Stainless_rev.htm
a = -1.4087;
b = 1.3982;
c = 0.2543;
d = -0.6260;
e = 0.2334;
f = 0.4256;
g = -0.4658;
h = 0.1650; 
i = -0.0199;

l = log10(T);

y = a + b * l + c * l.^2 + d * l.^3 + e * l.^4 + ...
        f * l.^5 + g * l.^6 + h * l.^7 + i* l.^8; 

eval = (10.^y);  

end
