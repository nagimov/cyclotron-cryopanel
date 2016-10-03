function [T_a, T_b] = RN_03
    clc; clear all;
    close all;
    tic

    % INTEGRATOR DATA
    t_delta = 0.1;  % time step, s
    t = 20;  % number of time steps, -
    HX_slices = 10;  % number of slices, -

    % HX DATA
    m = 1;  % mass flow, kg/s
    M = 1;  % mass of streams content within a cell, kg
    HX_UA = 10000;  % HX coefficient, W/K

    % INITIAL DATA
    fluid_a = 'helium';  % stream A fluid name
    p_a_in = 101325;  % inlet pressure of stream A, Pa
    T_a_in = 100;  % inlet temperature of stream A, K
    h_a_in = rp_htp(T_a_in, p_a_in, fluid_a);  % inlet enthalpy of stream A, J/kg
    fluid_b = 'nitrogen';  % stream B fluid name
    p_b_in = 101325;  % inlet pressure of stream B, Pa
    T_b_in = 200;  % inlet temperature of stream B, K
    h_b_in = rp_htp(T_b_in, p_b_in, fluid_b);  % inlet enthalpy of stream B, J/kg

    % INITIAL CONDITIONS
    h_a_0 = h_a_in*ones(1, HX_slices);
    h_b_0 = h_b_in*ones(1, HX_slices);
    p_a_0 = p_a_in*ones(1, HX_slices);
    p_b_0 = p_b_in*ones(1, HX_slices);
	h_sol(1, :) = [h_a_0    h_b_0];

	% SOLVER
	for j = 2 : t + 1
		[h_sol(j, :), fval, exitflag] = heateq(j);
        toc
        disp(['Time step completed ' num2str(j-1)...
        ' Time ' num2str(toc/60) ' min '...
        ' Fval ' num2str(sum(fval)) ...
        ' Exit Flag ' num2str(exitflag)])
    end
    
    %CONVERT TO T AND CALCULATE Q   
    T_a = zeros(t + 1, HX_slices);
    T_b = zeros(t + 1, HX_slices);

    for i = 1 : t + 1
        for j =  1 : HX_slices
            T_a(i,j) = rp_thp(h_sol(i,j), p_a_in, fluid_a);
            T_b(i,j) = rp_thp(h_sol(i,j + HX_slices), p_b_in, fluid_b);
            T_a_delta(i,j) = T_b(i,j) - T_a(i,j);
			T_b_delta(i,j) = -T_a_delta(i,j);
			Q_cond_a(i,j) = HX_UA / HX_slices * T_a_delta(i,j);
            Q_cond_b(i,j) = HX_UA / HX_slices * T_b_delta(i,j);
        end
    end     
   
    % HE PLOT 
    plot(1:HX_slices, T_a)
    xlabel('Slices')
    ylabel('Temperature')
    title('He')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_N2' num2str(HX_slices)],'-dpng','-r0')
    
    % N2 PLOT
    figure
    plot(1:HX_slices, T_b)
    xlabel('Slices')
    ylabel('Temperature')
    title('N_2')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_HE' num2str(HX_slices)],'-dpng','-r0')
    
    % OVERALL PLOT 
    figure 
    hold on 
    plot(1:HX_slices, T_a,'r')
    plot(1:HX_slices, T_b,'b')
    xlabel('Slices')
    ylabel('Temperature')
    title('Overall')
    axis([1 HX_slices T_a_in T_b_in])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_overall' num2str(HX_slices)],'-dpng','-r0')
    
    % OVERALL Q PLOT 
    figure 
    hold on 
    plot(1:HX_slices, Q_cond_a,'r')
    plot(1:HX_slices, Q_cond_b,'b')
    xlabel('Slices')
    ylabel('Heat (Q)')
    axis([1 HX_slices -1e5 1e5])
    title('Overall')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_Q' num2str(HX_slices)],'-dpng','-r0')    

	function [x, fval, exitflag] = heateq(j)  % time iteration
		options = optimset('TolX', 1e-7, 'TolFun', 1e-7, ...
		    			   'MaxFunEvals', 1e7, 'MaxIter', 1e7, ...
		    			   'Display', 'iter');
        h_a_prev = h_sol(j - 1, 1 : HX_slices);
        h_b_prev = h_sol(j - 1, HX_slices + 1 : 2 * HX_slices);
        guess = [h_a_prev h_b_prev];
		[x, fval, exitflag] = fsolve(@eqgen, guess, options);
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
            T_a = zeros(1, HX_slices);
			T_b = zeros(1, HX_slices);
			T_a_delta = zeros(1, HX_slices);
			T_b_delta = zeros(1, HX_slices);

			h_a = H(1 : HX_slices);
			h_b = H(HX_slices + 1 : 2 * HX_slices);
            p_a = p_a_0;
            p_b = p_b_0;

			% energy equations
			for i = 1 : HX_slices  % slicing loop

				T_a(i) = rp_thp(h_a(i), p_a(i), fluid_a);
				T_b(i) = rp_thp(h_b(i), p_b(i), fluid_b);

				T_a_delta(i) = T_b(i) - T_a(i);
				T_b_delta(i) = -T_a_delta(i);
				
				Q_cond_a(i) = HX_UA / HX_slices * T_a_delta(i);
                u_a(i) = rp_uhp(h_a(i), p_a(i), fluid_a);
				u_a_prev(i) = rp_uhp(h_a_prev(i), p_a(i), fluid_a);

				Q_cond_b(i) = HX_UA / HX_slices * T_b_delta(i);
				u_b(i) = rp_uhp(h_b(i), p_b(i), fluid_b);
				u_b_prev(i) = rp_uhp(h_b_prev(i), p_b(i), fluid_b);
                
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
end  % RN_03
