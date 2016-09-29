function RN_01
    clc; clear all;
    close all;
    tic
    global t_delta t HX_slices h_sol

    % INTEGRATOR DATA
    t_delta = 0.1;  % time step, s
    t = 5;  % number of time steps, -
    HX_slices = 5;  % number of slices, -

    % HX DATA
    m = 1;  % mass flow, kg/s
    M = 1;  % mass of streams content within a cell, kg
    HX_UA = 10000;  % HX coefficient, W/K

    % INITIAL DATA
    fluid_a = 'helium';  % stream A fluid name
    p_a_in = 101325;  % inlet pressure of stream A, Pa
    T_a_in = 100;  % inlet temperature of stream A, K
    h_a_in = refpropm('H', 'T', T_a_in, 'P', p_a_in, fluid_a);  % inlet enthalpy of stream A, J/kg
    u_a_in = refpropm('U', 'H', h_a_in, 'P', p_a_in, fluid_a);  % inlet internal energy of stream A, J/kg
    fluid_b = 'nitrogen';  % stream B fluid name
    p_b_in = 101325;  % inlet pressure of stream B, Pa
    T_b_in = 200;  % inlet temperature of stream B, K
    h_b_in = refpropm('H', 'T', T_b_in, 'P', p_b_in, fluid_b);  % inlet enthalpy of stream B, J/kg
    u_b_in = refpropm('U', 'H', h_b_in, 'P', p_b_in, fluid_b);  % inlet internal energy of stream B, J/kg

    % INITIAL CONDITIONS
    for i = 1 : HX_slices
        p_a_0(i) = p_a_in;  % Pa
        T_a_0(i) = T_a_in;  % K
        h_a_0(i) = refpropm('H', 'T', T_a_in, 'P', p_a_in, fluid_a);  % J/kg
        p_b_0(i) = p_b_in;  % Pa
        T_b_0(i) = T_b_in;  % K
        h_b_0(i) = refpropm('H', 'T', T_b_in, 'P', p_b_in, fluid_b);  % J/kg
	end
	h_sol(1, :) = [h_a_0    h_b_0];

	% SOLVER
	for j = 2 : t + 1
		[h_sol(j, :), fval(j, :), exitflag(j, :)] = heateq(j);
		toc
	end

	hold on
    plot(1:HX_slices, h_sol(:, 1 : HX_slices))
    xlabel('Slices')
    ylabel('Enthalpy')
    figure
    plot(1:HX_slices, h_sol(:, (HX_slices + 1 : 2 * HX_slices)))
    xlabel('Slices')
    ylabel('Enthalpy')

	function [x, fval, exitflag] = heateq(j)  % time iteration
		options = optimset('TolX', 1e-7, 'TolFun', 1e-7, ...
		    			   'MaxFunEvals', 1e7, 'MaxIter', 1e7, ...
		    			   'Display', 'iter');
		h_sol_prev = h_sol(j - 1, :);
        h_a_prev = h_sol_prev(1 : HX_slices);
        h_b_prev = h_sol_prev(HX_slices + 1 : 2 * HX_slices);
		[x, fval, exitflag] = fsolve(@eqgen, h_sol_prev, options);
		%[x, fval, exitflag] = fminsearch(@eqgen, h_sol_prev, options);
		% compute difference between the equation and zero
	    function F = eqgen(H)
	    	F_a = zeros(1, HX_slices);
	    	F_b = zeros(1, HX_slices);
	    	u_b = zeros(1, HX_slices);
	    	u_a = zeros(1, HX_slices);
	    	u_b_prev = zeros(1, HX_slices);
	    	u_a_prev = zeros(1, HX_slices);
			h_a_delta = zeros(1, HX_slices);
			h_b_delta = zeros(1, HX_slices);
			Q_cond_a = zeros(1, HX_slices);
			Q_cond_b = zeros(1, HX_slices);
			T_a_delta = zeros(1, HX_slices);
			T_b_delta = zeros(1, HX_slices);

			h_a = H(1 : HX_slices);
			h_b = H(HX_slices + 1 : 2 * HX_slices);

			% energy equations
			for i = 1 : HX_slices  % slicing loop
				% p_a_0
				% p_b_0
				T_a(i) = rp_thp(h_a(i), p_a_0(i), fluid_a);
				T_b(i) = rp_thp(h_b(i), p_b_0(i), fluid_b);

				T_a_delta(i) = T_b(i) - T_a(i);
				T_b_delta(i) = -T_a_delta(i);
				
				Q_cond_a(i) = HX_UA / HX_slices * T_a_delta(i);
				u_a(i) = rp_uhp(h_a(i), p_a_0(i), fluid_a);
				u_a_prev(i) = rp_uhp(h_a_prev(i), p_a_0(i), fluid_a);

				Q_cond_b(i) = HX_UA / HX_slices * T_b_delta(i);
				u_b(i) = rp_uhp(h_b(i), p_b_0(i), fluid_b);
				u_b_prev(i) = rp_uhp(h_b_prev(i), p_b_0(i), fluid_b);
			end
    		
    		% enthalpy deltas
    		for i = 1 : HX_slices
    			if i == 1
					h_a_delta(i) = h_a_in - h_a(1);
					h_b_delta(HX_slices) = h_b_in - h_b(HX_slices);
	            else
    				h_a_delta(i) = h_a(i - 1) - h_a(i);
    				h_b_delta(HX_slices + 1 - i) = h_b(HX_slices + 2 - i) - h_b(HX_slices + 1 - i);
	            end
			end

			for i = 1 : HX_slices
				F_a(i) = -(u_a(i) - u_a_prev(i)) * M / t_delta + (m * h_a_delta(i) + Q_cond_a(i));
				F_b(i) = -(u_b(i) - u_b_prev(i)) * M / t_delta + (m * h_b_delta(i) + Q_cond_b(i));
			end			
			% final exit vector
			F = [F_a    F_b];

		end
	end
end  % RN_01
