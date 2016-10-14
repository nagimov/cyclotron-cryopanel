function eval = cp(T)  % J/kg*K

% http://cryogenics.nist.gov/MPropsMAY/304Stainless/304Stainless_rev.htm
a = 22.0061;
b = -127.5528;
c = 303.647;
d = -381.0098;
e = 274.0328;
f = -112.9212;
g = 24.7593;
h = -2.239153; 

l = log10(T);

y = a + b * l + c * l.^2 + d * l.^3 + ...
    e * l.^4 + f * l.^5 + g * l.^6 + h * l.^7; 

eval = (10.^y);  

end
