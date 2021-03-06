function eval = propsc_htp(T, p, fluid, lib)
    [~, length] = size(T);
    buffer_size = 1000;
    ierr = 0;
    b = (1:1:buffer_size);
    herr = char(b);
    
    %creating input pointers
    input1Ptr = libpointer('doublePtr',p);
    input2Ptr = libpointer('doublePtr',T);
    
    %Selecting backend and fluid
    backend = 'BICUBIC&HEOS';
    [handle, sh] = calllib('coolprop','AbstractState_factory',backend,fluid,ierr,herr,buffer_size);
    [input_pair, sip] = calllib('coolprop','get_input_pair_index','PT_INPUTS');

    outputs=zeros(5,1);
    %Choosing parameters to compute
    [outputs(1,1), so1] = calllib('coolprop','get_param_index','Dmolar');
    [outputs(2,1), so2] = calllib('coolprop','get_param_index','Dmolar');
    [outputs(3,1), so3] = calllib('coolprop','get_param_index','Hmass');
    [outputs(4,1), so4] = calllib('coolprop','get_param_index','Smolar');
    [outputs(5,1), so5] = calllib('coolprop','get_param_index','Dmolar');
    
    %Creating ouput pointers
    out1Ptr = libpointer('doublePtr',zeros(length,1));
    out2Ptr = libpointer('doublePtr',zeros(length,1));
    out3Ptr = libpointer('doublePtr',zeros(length,1));
    out4Ptr = libpointer('doublePtr',zeros(length,1));
    out5Ptr = libpointer('doublePtr',zeros(length,1));
    
    calllib('coolprop','AbstractState_update_and_5_out',handle,input_pair,input1Ptr,input2Ptr,length,outputs,out1Ptr,out2Ptr,out3Ptr,out4Ptr,out5Ptr,ierr,herr,buffer_size);

    %Saving computed values to array
    out1=get(out1Ptr,'Value');
    out2=get(out2Ptr,'Value');
    out3=get(out3Ptr,'Value');
    out4=get(out4Ptr,'Value');
    out5=get(out5Ptr,'Value');

    eval=out3';
end
