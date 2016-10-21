function [T_a, T_b] = RN_03n % Most basic model 
    clc; clear all;
    close all;
    tic

    % COOLPROP
    CP_dump = 20; % number of time steps before we dump & re-load the library 
    path_to_lib = 'D:\CoolProp_wrapper_fast'; %specify path to coolprop shared library
    path_to_include= 'D:\CoolProp_wrapper_fast'; %specify path to coolprop's include folder
    libname = 'libCoolProp'; % OSX and linux
        if ispc
            libname = 'CoolProp';
        end
    addpath(path_to_lib)
    addpath(path_to_include)
    loadcoolprop; 
       
    % Loading shared library
    function loadcoolprop 
        if ~libisloaded('coolprop') %checking whether library is already loaded
            loadlibrary(libname,'CoolPropLib.h','includepath',...
                path_to_include,'alias','coolprop'); % loading library with alias coolprop
        end
    end

    % INTEGRATOR DATA
    t_delta = 0.1;  % time step, s
    t = 20;  % number of time steps, -
    HX_slices = 20;  % number of slices, -

    % HX DATA
    m = 1;  % mass flow, kg/s
    M = 1;  % mass of streams content within a cell, kg
    HX_UA = 2500;  % HX coefficient, W/K

    % INITIAL DATA
    fluid_a = 'helium';  % stream A fluid name
    p_a_in = 101325;  % inlet pressure of stream A, Pa
    q_a_in = .1;  % inlet temperature of stream A, K
    h_a_in = prop_hqp(q_a_in, p_a_in, fluid_a, 'CP');  % inlet enthalpy of stream A, J/kg
    fluid_b = 'nitrogen';  % stream B fluid name
    p_b_in = 101325;  % inlet pressure of stream B, Pa
    T_b_in = 100;  % inlet temperature of stream B, K
    h_b_in = prop_htp(T_b_in, p_b_in, fluid_b, 'CP');  % inlet enthalpy of stream B, J/kg
    T_a_in = 4; % Used only as a lower T limit for y-axis on plots


    % INITIAL CONDITIONS
    h_a_0 = h_a_in * ones(HX_slices, 1);
    h_b_0 = h_b_in * ones(HX_slices, 1);
    p_a_0 = p_a_in * ones(HX_slices, 1);
    p_b_0 = p_b_in * ones(HX_slices, 1);
	h_a_sol(:, 1) = h_a_0;
    h_b_sol(:, 1) = h_b_0;

	% SOLVER
	for j = 2 : t + 1
		[sol, fval, exitflag] = heateq(j);
        disp(['Time step completed ' num2str(j-1)...
        ' Time ' num2str(toc/60) ' min '...
        ' Fval ' num2str(sum(fval)) ...
        ' Exit Flag ' num2str(exitflag)])
        h_a_sol(:, j) = sol(:, 1);
        h_b_sol(:, j) = sol(:, 2);
        
        % dump and re-load coolprop
        if rem(j, CP_dump) == 0
        unloadlibrary 'coolprop'
        loadcoolprop; 
        end 
    end
    
    % CONVERT TO T  
    T_a = zeros(HX_slices, t + 1);
    T_b = zeros(HX_slices, t + 1);
    
    for j =  1 : t + 1
        T_a(:, j) = props_thp(h_a_sol(:, j), p_a_0, fluid_a, 'CP');
        T_b(:, j) = props_thp(h_b_sol(:, j), p_b_0, fluid_b, 'CP');
    end
   
    % HE PLOT 
    plot(1:HX_slices, T_a, 'r')
    xlabel('Slices')
    ylabel('Temperature')
    title('He')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_N2' num2str(HX_slices)],'-dpng','-r0')
    
    % N2 PLOT
    figure
    plot(1:HX_slices, T_b, 'b')
    xlabel('Slices')
    ylabel('Temperature')
    title('N_2')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_HE' num2str(HX_slices)],'-dpng','-r0')
    
    % OVERALL PLOT 
    figure 
    hold on 
    plot(1:HX_slices, T_a, 'r')
    plot(1:HX_slices, T_b, 'b')
    xlabel('Slices')
    ylabel('Temperature')
    title('Overall')
    axis([1 HX_slices T_a_in T_b_in])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_overall' num2str(HX_slices)],'-dpng','-r0')

    %*********************************************************************
    % FUNCTION THAT CALCULATES Q
    function Q_cond = Q_cond(h_a, p_a, h_b, p_b)
        T_a = propsc_thp(h_a, p_a', fluid_a, 'CP');
        T_b = propsc_thp(h_b, p_b', fluid_b, 'CP');
        T_delta = T_b - T_a;
        Q_cond = HX_UA / HX_slices * T_delta;
    end
    
    % FUNCTION THAT COMPUTES du_dt
    function dudt = dudt(h, h_prev, p, fluid)
		dudt = propsc_uhp(h, p', fluid, 'CP') - propsc_uhp(h_prev, p', fluid, 'CP');
    end
    
    %*********************************************************************
	function [x, fval, exitflag] = heateq(j)  % time iteration
		options = optimset('TolX', 1e-7, 'TolFun', 1e-7, ...
		    			   'MaxFunEvals', 1e7, 'MaxIter', 1e7, ...
		    			   'Display', 'iter');
        h_a_prev = h_a_sol(:, j - 1);
        h_b_prev = h_b_sol(:, j - 1);
        guess = [h_a_prev h_b_prev];
		[x, fval, exitflag] = fsolve(@eqgen, guess, options);

        % compute difference between the equation and zero
	    function eps = eqgen(sol)

            % pre-allocate 
			h_a_delta = zeros(HX_slices, 1);
			h_b_delta = zeros(HX_slices, 1);

            % split solution vector 
			h_a = sol(:, 1);
			h_b = sol(:, 2);
            p_a = p_a_0;
            p_b = p_b_0;

			% energy equations
            dudt_a = dudt(h_a, h_a_prev, p_a, fluid_a);
            dudt_b = dudt(h_b, h_b_prev, p_b, fluid_b);
			Q_cond_a = Q_cond(h_a, p_a, h_b, p_b);
            Q_cond_b = -Q_cond_a; 

    		% enthalpy deltas
    		for i = 1 : HX_slices
    			if i == 1
					h_a_delta(1) = h_a_in - h_a(1);
					h_b_delta(HX_slices) = h_b_in - h_b(HX_slices);
	            else
    				h_a_delta(i) = h_a(i - 1) - h_a(i);
    				h_b_delta(HX_slices + 1 - i) = h_b(HX_slices + 2 - i) - h_b(HX_slices + 1 - i);
	            end
            end
            
            % final equations 
			eps_a = -dudt_a * M / t_delta + (m * h_a_delta + Q_cond_a);
			eps_b = -dudt_b * M / t_delta + (m * h_b_delta + Q_cond_b);
                
			% final exit vector
			eps = [eps_a  eps_b];
		end
	end
end  % RN_03n 
