function fval =  refpropwrap(Prop,H,P,Name) 

try
    fval = refpropm(Prop,'H',H,'P',P,Name);  
catch ME 
    %Pick out numbers from error message
    val = regexp(ME.message, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?','match');
    
    if str2double(val{1}) == 249 %Code 249 means the error is on H limits  
    h_min = str2double(val{3}); %Get h_min value 
    h_max = str2double(val{4}); %Get h_max value 
    
    if H < h_min
      fval = refpropwrap(Prop,h_min,P,Name);
      display('using h_min value!')
    else
      fval = refpropwrap(Prop,h_max,P,Name);
      display('using h_max value!')
    end
    
    elseif str2double(val{1}) == 4 %Code 4 means the error is on the upper P limit
      p_max = str2double(val{4}); %Get p_max value 
      fval = refpropwrap(Prop,H,p_max,Name);
      display('using p_max value!')
    else   
    end
end 
