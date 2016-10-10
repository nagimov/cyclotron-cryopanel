function f =  prop_hqp(q,p,Name,lib) 

if strcmp(lib,'RP') == 1
    p = p/1000;
    f = refpropm('H','P',p,'Q',q,Name);  
elseif strcmp(lib,'CP') == 1
    f = CoolProp.PropsSI('H','P',p,'Q',q,Name);
end

end 
