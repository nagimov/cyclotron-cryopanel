function f =  rp_hpq(p,q,Name) 

p = p/1000;

f = refpropm('H','P',p,'Q',q,Name);  

end 
