function [T_a_sol, T_b_sol, T_w_sol] = RN_04 %With wall 
    addpath('H:\coolprop_lib')
    clc; clear all;
    close all;
    tic

    % INTEGRATOR DATA
    t_delta = 0.1;  % time step, s
    t = 5;  % number of time steps, -
    HX_slices = 10;  % number of slices, -
    lib = 'CP'; 

    % HX DATA
    m = 1;  % mass flow, kg/s
    M = 1;  % mass of streams content within a cell, kg
    M_w = 1; % mass of wall, kg 
    HX_UA_data = {'nitrogen', 10000; ...
                    'helium', 10000}; % HX coefficient, W/K

    % INITIAL DATA
    fluid_a = 'helium';  % stream A fluid name
    p_a_in = 101325;  % inlet pressure of stream A, Pa
    T_a_in = 100;  % inlet temperature of stream A, K
    h_a_in = prop_htp(T_a_in, p_a_in, fluid_a, lib);  % inlet enthalpy of stream A, J/kg
    fluid_b = 'nitrogen';  % stream B fluid name
    p_b_in = 101325;  % inlet pressure of stream B, Pa
    T_b_in = 200;  % inlet temperature of stream B, K
    h_b_in = prop_htp(T_b_in, p_b_in, fluid_b, lib);  % inlet enthalpy of stream B, J/kg
    T_w_in = 150; % intitial wall temperature, K

    % INITIAL CONDITIONS
    h_a_0 = h_a_in * ones(1, HX_slices);
    h_b_0 = h_b_in * ones(1, HX_slices);
    p_a_0 = p_a_in * ones(1, HX_slices);
    p_b_0 = p_b_in * ones(1, HX_slices);
    T_w_0 = T_w_in * ones(1, HX_slices); 
	h_a_sol(1, :) = h_a_0;
    h_b_sol(1, :) = h_b_0; 
    T_w_sol(1, :) = T_w_0;

	% SOLVER
	for j = 2 : t + 1
		[sol, fval, exitflag] = heateq(j);
		disp(['Time step completed ' num2str(j-1)...
            ' Time ' num2str(toc/60) ' min '...
            ' Fval ' num2str(sum(fval)) ...
            ' Exit Flag ' num2str(exitflag)])
        h_a_sol(j, :) = sol(1 : 1 * HX_slices );
        h_b_sol(j, :) = sol(1 * HX_slices + 1 : 2 * HX_slices );
        T_w_sol(j, :) = sol(2 * HX_slices + 1 : 3 * HX_slices );
    end
    
    % FUNCTION THAT PULLS HX_UA DATA
    function value = HX_UA(fluid)        
        f = strcmp(fluid, HX_UA_data);   
        value = HX_UA_data{f, 2};
    end

    % FUNCTION THAT CALCULATES Q
    function Q_cond = Q_cond(h, p, T_w, fluid)
        T = prop_thp(h, p, fluid, lib);
        T_delta = T_w - T;
        Q_cond = HX_UA(fluid) / HX_slices * T_delta;
    end 

    % FUNCTION TAHT COMPUTES delta_u
    function delta_u = u_delta(h, h_prev, p, fluid)
		delta_u = prop_uhp(h, p, fluid, lib) - prop_uhp(h_prev, p, fluid, lib); 
    end
        
    % CONVERT TO T    
    T_a_sol = zeros(t + 1, HX_slices);
    T_b_sol = zeros(t + 1, HX_slices);

    for i = 1 : t + 1
        for j =  1 : HX_slices
            T_a_sol(i,j) = prop_thp(h_a_sol(i,j), p_a_in, fluid_a, lib);
            T_b_sol(i,j) = prop_thp(h_b_sol(i,j), p_b_in, fluid_b, lib);
        end
    end
     
    % HE PLOT 
    plot(1:HX_slices, T_a_sol)
    xlabel('Slices')
    ylabel('Temperature')
    title('He')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_N2' num2str(HX_slices)],'-dpng','-r0')
    
    % N2 PLOT
    figure
    plot(1:HX_slices, T_b_sol)
    xlabel('Slices')
    ylabel('Temperature')
    title('N_2')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_HE' num2str(HX_slices)],'-dpng','-r0')
    
    % OVERALL PLOT
    figure 
    hold on 
    plot(1:HX_slices, T_a_sol, 'r')
    plot(1:HX_slices, T_b_sol, 'b')
    plot(1:HX_slices, T_w_sol, 'g')
    xlabel('Slices')
    ylabel('Temperature')
    title('Overall')
    axis([1 HX_slices T_a_in T_b_in])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_overall' num2str(HX_slices)],'-dpng','-r0')
    
    %*********************************************************************

	function [x, fval, exitflag] = heateq(j)  % time iteration
		options = optimset('TolX', 1e-7, 'TolFun', 1e-7, ...
		    			   'MaxFunEvals', 1e7, 'MaxIter', 1e7, ...
		    			   'Display', 'iter');
        h_a_prev = h_a_sol(j - 1, :);
        h_b_prev = h_b_sol(j - 1, :);
        T_w_prev = T_w_sol(j - 1, :); 
        sol_guess = [h_a_prev h_b_prev T_w_prev];
		[x, fval, exitflag] = fsolve(@eqgen, sol_guess, options);

        % compute difference between the equation and zero
	    function F = eqgen(sol)
            % pre-allocate 
	    	F_a = zeros(1, HX_slices);
	    	F_b = zeros(1, HX_slices);
   	    	F_w = zeros(1, HX_slices);
	    	u_a_delta = zeros(1, HX_slices);
	    	u_b_delta = zeros(1, HX_slices);
			h_a_delta = zeros(1, HX_slices);
			h_b_delta = zeros(1, HX_slices);
			Q_cond_aw = zeros(1, HX_slices);
			Q_cond_bw = zeros(1, HX_slices);
            T_w_delta = zeros(1, HX_slices); 
            cp_w = zeros(1, HX_slices); 
            
            % split solution vector 
			h_a = sol(1 : 1 * HX_slices);
			h_b = sol(1 * HX_slices + 1 : 2 * HX_slices);
            T_w = sol(2 * HX_slices + 1 : 3 * HX_slices);
            p_a = p_a_0;
            p_b = p_b_0;

			% energy equations
			for i = 1 : HX_slices  % slicing loop                
                Q_cond_aw(i) = Q_cond(h_a(i), p_a(i), T_w(i), fluid_a);
                Q_cond_bw(i) = Q_cond(h_b(i), p_b(i), T_w(i), fluid_b);
                u_a_delta(i) = u_delta(h_a(i), h_a_prev(i), p_a(i), fluid_a);
                u_b_delta(i) = u_delta(h_b(i), h_b_prev(i), p_b(i), fluid_b);
                T_w_delta(i) = T_w(i) - T_w_prev(i); 
                cp_w(i) = cp(T_w(i)); 
                        
            % enthalpy deltas  
    			if i == 1
					h_a_delta(1) = h_a_in - h_a(1);
					h_b_delta(HX_slices) = h_b_in - h_b(HX_slices);
	            else
    				h_a_delta(i) = h_a(i - 1) - h_a(i);
    				h_b_delta(HX_slices + 1 - i) = h_b(HX_slices + 2 - i) - h_b(HX_slices + 1 - i);
	            end
            end
            
            % final equations 
			for i = 1 : HX_slices
				F_a(i) = - u_a_delta(i) * M / t_delta + m * h_a_delta(i) + Q_cond_aw(i);
				F_b(i) = - u_b_delta(i) * M / t_delta + m * h_b_delta(i) + Q_cond_bw(i);
                F_w(i) = - T_w_delta(i) * M_w / t_delta + (Q_cond_aw(i) + Q_cond_bw(i)) / cp_w(i) ;
            end
            
			% combine final exit vector
			F = [F_a  F_b  F_w];
		end
	end
end  % RN_04
