addpath('D:\CoolProp_wrapper\')
clear all;
close all;

N = 1e4;
M = 10;

lib = 'CP';
fluid = 'helium';
p_test = linspace(101325, 2 * 101325, N);
T_test = linspace(200, 300, N);
times = zeros(1, M);

tic
for i = 1 : M
    h_test_1 = props_htp(T_test, p_test, fluid, lib);
    times(i) = toc;
end
disp(mean(times))


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


tic
for i = 1 : M
    h_test_2 = propsc_htp(T_test, p_test, fluid, lib);
    times(i) = toc;
end
disp(mean(times))

delta = sum(abs(h_test_1 - h_test_2)) / sum(abs(h_test_1))