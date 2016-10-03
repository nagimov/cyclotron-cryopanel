function [T_a_sol, T_b_sol, T_w_sol] = RN_05 %Time step adjustement
    clc; clear all;
    close all;
    tic

    % INTEGRATOR DATA
    t_delta = 0.1;  % time step, s
    t = 20;  % number of time steps, -
    HX_slices = 10;  % number of slices, -
    N = 2; % frequency of calcuing Q 

    % HX DATA
    m = 1;  % mass flow, kg/s
    M = 1;  % mass of streams content within a cell, kg
    M_w = 1; % mass of wall, kg 
    HX_UA_aw = 10000;  % HX coefficient, stream A to wall, W/K
    HX_UA_bw = 10000;  % HX coefficient, stream B to wall, W/K

    % INITIAL DATA
    fluid_a = 'helium';  % stream A fluid name
    p_a_in = 101325;  % inlet pressure of stream A, Pa
    T_a_in = 100;  % inlet temperature of stream A, K
    h_a_in = rp_htp(T_a_in, p_a_in, fluid_a);  % inlet enthalpy of stream A, J/kg
    fluid_b = 'nitrogen';  % stream B fluid name
    p_b_in = 101325;  % inlet pressure of stream B, Pa
    T_b_in = 200;  % inlet temperature of stream B, K
    h_b_in = rp_htp(T_b_in, p_b_in, fluid_b);  % inlet enthalpy of stream B, J/kg
    T_w_in = 150; 

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
    
    %CONVERT TO T    
    T_a_sol = zeros(t + 1, HX_slices);
    T_b_sol = zeros(t + 1, HX_slices);

    for i = 1 : t + 1
        for j =  1 : HX_slices
            T_a_sol(i,j) = rp_thp(h_a_sol(i,j), p_a_in, fluid_a);
            T_b_sol(i,j) = rp_thp(h_b_sol(i,j), p_b_in, fluid_b);
        end
    end  
    
    % HE PLOT 
    plot(1:HX_slices, T_a_sol, 'r')
    xlabel('Slices')
    ylabel('Temperature')
    title('He')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_N2' num2str(HX_slices)],'-dpng','-r0')
    
    % N2 PLOT
    figure
    plot(1:HX_slices, T_b_sol, 'b')
    xlabel('Slices')
    ylabel('Temperature')
    title('N_2')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_HE' num2str(HX_slices)],'-dpng','-r0')
    
    % Make Odd and Even vector for plots (so far works for only N=2!!)
    t_all = 1 : t + 1;
    t_x = 1 : N : t + 1; 
    t_y = setdiff(t_all, t_x);
    
    % OVERALL PLOT
    figure 
    hold on 
    plot(1:HX_slices, T_a_sol(t_x, :), 'r:', 'LineWidth', 2)
    plot(1:HX_slices, T_a_sol(t_y, :)', 'r', 'LineWidth', 2)
    plot(1:HX_slices, T_b_sol(t_x, :), 'b:', 'LineWidth', 2)
    plot(1:HX_slices, T_b_sol(t_y, :)', 'b', 'LineWidth', 2)
    plot(1:HX_slices, T_w_sol(t_x, :), 'g:', 'LineWidth', 2)
    plot(1:HX_slices, T_w_sol(t_y, :)', 'g', 'LineWidth', 2)
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
        
        function Q_cond = Q_cond_calc(h, p, T_w, fluid)
            T = rp_thp(h, p, fluid);
            T_delta = T_w - T;
            Q_cond = HX_UA_aw / HX_slices * T_delta;
        end    
        
        % compute difference between the equation and zero
	    function F = eqgen(sol)
            % pre-allocate 
	    	F_a = zeros(1, HX_slices);
	    	F_b = zeros(1, HX_slices);
   	    	F_w = zeros(1, HX_slices);
	    	u_b = zeros(1, HX_slices);
	    	u_a = zeros(1, HX_slices);
	    	u_b_prev = zeros(1, HX_slices);
	    	u_a_prev = zeros(1, HX_slices);
			h_a_delta = zeros(1, HX_slices);
			h_b_delta = zeros(1, HX_slices);
			Q_cond_aw = zeros(1, HX_slices);
			Q_cond_bw = zeros(1, HX_slices);
            Q_cond_aw_prev = zeros(1, HX_slices);
			Q_cond_bw_prev = zeros(1, HX_slices);
            cp_w = zeros(1, HX_slices); 
            
            % split solution vector 
			h_a = sol(1 : 1 * HX_slices);
			h_b = sol(1 * HX_slices + 1 : 2 * HX_slices);
            T_w = sol(2 * HX_slices + 1 : 3 * HX_slices);
            p_a = p_a_0;
            p_b = p_b_0;

			% energy equations
			for i = 1 : HX_slices  % slicing loop
                
                if rem(j, N) == 0
                    Q_cond_aw(i) = Q_cond_calc(h_a(i),p_a(i),T_w(i),fluid_a);
                    Q_cond_bw(i) = Q_cond_calc(h_b(i),p_b(i),T_w(i),fluid_b);
                    Q_cond_aw_prev(i) = Q_cond_aw(i);
                    Q_cond_bw_prev(i) = Q_cond_bw(i);
                else
                    Q_cond_aw(i) = Q_cond_aw_prev(i);
                    Q_cond_bw(i) = Q_cond_bw_prev(i);
                end     

                cp_w(i) = cp(T_w(i));                 
				u_a(i) = rp_uhp(h_a(i), p_a(i), fluid_a);
				u_a_prev(i) = rp_uhp(h_a_prev(i), p_a(i), fluid_a);
				u_b(i) = rp_uhp(h_b(i), p_b(i), fluid_b);   
     			u_b_prev(i) = rp_uhp(h_b_prev(i), p_b(i), fluid_b);             
           
            % enthalpy deltas  
    			if i == 1
					h_a_delta(i) = h_a_in - h_a(1);
					h_b_delta(HX_slices) = h_b_in - h_b(HX_slices);
	            else
    				h_a_delta(i) = h_a(i - 1) - h_a(i);
    				h_b_delta(HX_slices + 1 - i) = h_b(HX_slices + 2 - i) - h_b(HX_slices + 1 - i);
	            end
            end
            
            % final equations 
			for i = 1 : HX_slices
				F_a(i) = -(u_a(i) - u_a_prev(i)) * M / t_delta + (m * h_a_delta(i) + Q_cond_aw(i));
				F_b(i) = -(u_b(i) - u_b_prev(i)) * M / t_delta + (m * h_b_delta(i) + Q_cond_bw(i));
                F_w(i) = -(T_w(i) - T_w_prev(i)) * M_w / t_delta + (Q_cond_aw(i) + Q_cond_bw(i)) / cp_w(i) ;
            end
            
			% combine final exit vector
			F = [F_a  F_b  F_w];
		end
	end
end  % RN_05
