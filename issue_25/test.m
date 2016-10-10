%% Calling low-level interface from MATLAB through shared library (DLL)
%%
%% Originally developed by Edo Macchi, modified by Ian Bell
%% 
%% December 2015

path_to_lib = 'D:\CoolProp_wrapper_fast'; %specify path to coolprop shared library
path_to_include= 'D:\CoolProp_wrapper_fast'; %specify path to coolprop's include folder

% Loading shared library
if ~libisloaded('coolprop') %checking whether library is already loaded
    addpath(path_to_lib)
    addpath(path_to_include)
    libname = 'libCoolProp' % OSX and linux
    if ispc
        libname = 'CoolProp'
    end
    loadlibrary(libname,'CoolPropLib.h','includepath',path_to_include,'alias','coolprop'); % loading library with alias coolprop
    disp('loaded CoolProp shared library.')
    disp('loaded these functions: ')
    libfunctions coolprop
end

buffer_size = 1000;
ierr = 0;
b = (1:1:buffer_size);
herr= char(b);

%Selecting backend and fluid
backend = 'BICUBIC&HEOS';
fluid = 'Helium';
[handle, sh] = calllib('coolprop','AbstractState_factory',backend,fluid,ierr,herr,buffer_size);

length=100000;
input2 = linspace(200.0,300.0,length)';
input1 = linspace(101325,2*101325,length);
input1Ptr = libpointer('doublePtr',input1);
input2Ptr = libpointer('doublePtr',input2);
[input_pair,sip] = calllib('coolprop','get_input_pair_index','PT_INPUTS');

outputs=zeros(5,1);
%Choosing parameters to compute
[outputs(1,1), so1] = calllib('coolprop','get_param_index','T');
[outputs(2,1), so2] = calllib('coolprop','get_param_index','P');
[outputs(3,1), so3] = calllib('coolprop','get_param_index','Dmolar');
[outputs(4,1), so4] = calllib('coolprop','get_param_index','Hmolar');
[outputs(5,1), so5] = calllib('coolprop','get_param_index','Smolar');

%Creating ouput pointers
out1Ptr = libpointer('doublePtr',zeros(length,1));
out2Ptr = libpointer('doublePtr',zeros(length,1));
out3Ptr = libpointer('doublePtr',zeros(length,1));
out4Ptr = libpointer('doublePtr',zeros(length,1));
out5Ptr = libpointer('doublePtr',zeros(length,1));

tic;

calllib('coolprop','AbstractState_update_and_5_out',handle,input_pair,input1Ptr,input2Ptr,length,outputs,out1Ptr,out2Ptr,out3Ptr,out4Ptr,out5Ptr,ierr,herr,buffer_size);

dt=toc;

%Saving computed values to array
out1=get(out1Ptr,'Value');
out2=get(out2Ptr,'Value');
out3=get(out3Ptr,'Value');
out4=get(out4Ptr,'Value');
out5=get(out5Ptr,'Value')

fprintf('%d us/call\n', (dt*1e6/length));