function eval =  props_uhp(h,p,Name,lib) 

% LOOP 
if size(h) == size(p)   %Check to ensure input matrices are the same size
    
    [N, M] = size(h);
    eval = zeros(N, M); 
    for i = 1 : N
        for j = 1 : M
            eval(i,j) = prop_uhp(h(i,j), p(i,j), Name, lib);
        end
    end     
        
else
    error('Size of input matrices for u and p must be the same')
end 


% ACTUAL PROP FUNCTION 
function f = prop_uhp(h,p,Name,lib)

if strcmp(lib,'RP') == 1
    try
        p = p/1000;
        f = refpropm('U','H',h,'P',p,Name);  
    catch ME 
        val = regexp(ME.message, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?','match');
    
    if str2double(val{1}) == 249 %Code 249 means the error is on H limits  
        h_min = str2double(val{3}); %Get h_min value 
        h_max = str2double(val{4}); %Get h_max value 
    
    if h < h_min
      f = rp_uhp(h_min,p,Name);
      display('using h_min value!')
    else
      f = rp_uhp(h_max,p,Name);
      display('using h_max value!')
    end
    
    elseif str2double(val{1}) == 4 && p > 0 %Code 4 means the error is on the upper P limit
      p_max = str2double(val{4}); %Get p_max value 
      f = rp_uhp(h,p_max,Name);
      display('using p_max value!')
    elseif p < 0
      p_min = 0;  
      f = rp_uhp(h,p_min,Name);
      display('using p_min value!')  
    else   
      disp(ME.message)
      f = 'unknown error, see above msg'; 
    end
    end
    
elseif strcmp(lib,'CP') == 1
    try 
        f = CoolProp.PropsSI('U','H',h,'P',p,Name);
    catch ME
    end
end 
