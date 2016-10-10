function eval =  props_htp(T,p,Name,lib) 

% LOOP 
if size(T) == size(p)   %Check to ensure input matrices are the same size
    
    [N, M] = size(T);
    eval = zeros(N, M); 
    for i = 1 : N
        for j = 1 : M
            eval(i,j) = prop_htp(T(i,j), p(i,j), Name, lib);
        end
    end     
        
else
    error('Size of input matrices for T and p must be the same')
end 
end

% ACTUAL PROP FUNCTION 
function f = prop_htp(T,p,Name,lib)
    if strcmp(lib,'RP') == 1
        p = p/1000;
        f = refpropm('H','T',T,'P',p,Name);  
    elseif strcmp(lib,'CP') == 1
        f = CoolProp.PropsSI('H','T',T,'P',p,Name);  
    end
end
    
