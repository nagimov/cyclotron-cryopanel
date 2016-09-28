function f =  rp_thp(h,p,Name) 

p = p/1000;

try
    f = refpropm('T','H',h,'P',p,Name);  
catch ME 
    %Pick out numbers from error message
    val = regexp(ME.message, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?','match');
    
    if str2double(val{1}) == 249 %Code 249 means the error is on H limits  
    h_min = str2double(val{3}); %Get h_min value 
    h_max = str2double(val{4}); %Get h_max value 
    
    if h < h_min
      f = rp_thp(h_min,p,Name);
      display('using h_min value!')
    else
      f = rp_thp(h_max,p,Name);
      display('using h_max value!')
    end
    
    elseif str2double(val{1}) == 4 && p > 0 %Code 4 means the error is on the upper P limit
      p_max = str2double(val{4}); %Get p_max value 
      f = rp_thp(h,p_max,Name);
      display('using p_max value!')
    elseif p < 0
      p_min = 0;  
      f = rp_thp(h,p_min,Name);
      display('using p_min value!')  
    else   
      disp(ME.message)
      f = 'unknown error, see above msg'; 
    end
end 
