function [L,mu,rho,Pr,cp] = propsc_LmurhoPr_hp(h, p, fluid, ~) 
    
    [a, b] = size(h);
    length = max(a,b); % essentially determines if the input is 
                       % a column or row vector
                       
    buffer_size = 1000;
    ierr = 0;
    b = (1:1:buffer_size);
    herr = char(b);
    
    %creating input pointers
    input1Ptr = libpointer('doublePtr',h);
    input2Ptr = libpointer('doublePtr',p);
    
    %Selecting backend and fluid
    backend = 'BICUBIC&HEOS';
    [handle, ~] = calllib('coolprop','AbstractState_factory',backend,fluid,...
        ierr,herr,buffer_size);
    [input_pair, ~] = calllib('coolprop','get_input_pair_index','HmassP_INPUTS');

    outputs=zeros(5,1);
    %Choosing parameters to compute
    [outputs(1,1), ~] = calllib('coolprop','get_param_index','L');
    [outputs(2,1), ~] = calllib('coolprop','get_param_index','V');
    [outputs(3,1), ~] = calllib('coolprop','get_param_index','Dmass');
    [outputs(4,1), ~] = calllib('coolprop','get_param_index','Prandtl');
    [outputs(5,1), ~] = calllib('coolprop','get_param_index','Cpmass');
    
    %Creating ouput pointers
    out1Ptr = libpointer('doublePtr',zeros(length,1));
    out2Ptr = libpointer('doublePtr',zeros(length,1));
    out3Ptr = libpointer('doublePtr',zeros(length,1));
    out4Ptr = libpointer('doublePtr',zeros(length,1));
    out5Ptr = libpointer('doublePtr',zeros(length,1));
    
    calllib('coolprop','AbstractState_update_and_5_out',...
        handle,input_pair,input1Ptr,input2Ptr,length,outputs,...
        out1Ptr,out2Ptr,out3Ptr,out4Ptr,out5Ptr,...
        ierr,herr,buffer_size);

    %Getting output 
    L=get(out1Ptr,'Value'); % W/(m * K)
    mu=get(out2Ptr,'Value'); % Pa * s
    rho=get(out3Ptr,'Value'); % kg/m^3
    Pr=get(out4Ptr,'Value'); 
    cp=get(out5Ptr,'Value'); 

end
