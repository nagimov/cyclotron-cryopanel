function fval =  refpropwrap(Prop,T,P,Name) 

try
    fval = refpropm(Prop,'T',T,'P',P,Name); %If all limits are OK, just evaluate  
catch
    fval = 1e10; 
end 
