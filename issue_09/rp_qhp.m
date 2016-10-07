function f =  rp_qhp(h,p,Name) 

p = p/1000;

f = refpropm('Q','H',h,'P',p,Name);  

end 
