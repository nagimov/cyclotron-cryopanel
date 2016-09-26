function fval =  refpropwrap(Name,Prop,T,P) 

a = strcmp(Name,'nitrogen'); %Check to see if input matches nitrogen 
b = strcmp(Name,'helium'); %Check to see if input matches helium 

if a == 1 && T <= 63.151 %If a was indeed a positive comparison, then check limit for nitrogen 
    fval = 1e10;
elseif b == 1 && T <= 2.1768 %If b was indeed a positive comparison, then check limit for helium 
    fval = 1e10;
else 
    fval = refpropm(Prop,'T',T,'P',P,Name); %If all limits are OK, just evaluate  
end 



