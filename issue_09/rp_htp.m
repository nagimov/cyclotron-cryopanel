function f =  rp_htp(T,p,Name) 

p = p/1000;

f = refpropm('H','T',T,'P',p,Name);  

end 
