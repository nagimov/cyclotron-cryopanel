function [h_sol, T_sol] = RN_02
    clc; clear;
    close all;
    tic

    % INTEGRATOR DATA
    t_delta = 0.1;  % time step, s
    t = 5;  % number of time steps, -
    HX_slices = 5;  % number of slices, -

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
    h_a_in = refpropm('H', 'T', T_a_in, 'P', p_a_in/1e3, fluid_a);  % inlet enthalpy of stream A, J/kg
    fluid_b = 'nitrogen';  % stream B fluid name
    p_b_in = 101325;  % inlet pressure of stream B, Pa
    T_b_in = 200;  % inlet temperature of stream B, K
    h_b_in = refpropm('H', 'T', T_b_in, 'P', p_b_in/1e3, fluid_b);  % inlet enthalpy of stream B, J/kg
    T_w_in = 150; 

    % INITIAL CONDITIONS
    h_a_0 = h_a_in*ones(1, HX_slices);
    h_b_0 = h_b_in*ones(1, HX_slices);
    p_a_0 = p_a_in*ones(1, HX_slices);
    p_b_0 = p_b_in*ones(1, HX_slices);
    T_w_0 = T_w_in*ones(1, HX_slices); 
	h_sol(1, :) = [h_a_0    h_b_0];
    T_sol(1, :) = T_w_0;

	% SOLVER
	for j = 2 : t + 1
		[sol(j, :), fval, exitflag] = heateq(j);
		disp(['Time step completed ' num2str(j-1)...
            ' Time ' num2str(toc/60) ' min '...
            ' Fval ' num2str(sum(fval)) ...
            ' Exit Flag ' num2str(exitflag)])
        h_sol(j, :) = sol(j, 1 : 2 * HX_slices );
        T_sol(j, :) = sol(j, 2 * HX_slices + 1 : 3 * HX_slices );
    end

    %CONVERT TO T    
    for i = 1 : t + 1
        for j =  1 : HX_slices
            T(i,j) = refpropm('T','H',h_sol(i,j),'P',p_a_in/1e3,fluid_a);
        end
    end
    
    for i = 1 : t + 1
        for j = HX_slices + 1 : 2 * HX_slices 
            T(i,j) = refpropm('T','H',h_sol(i,j),'P',p_a_in/1e3,fluid_b);
        end
    end  
    
    % HE PLOT 
    plot(1:HX_slices, T(:, 1 : HX_slices))
    xlabel('Slices')
    ylabel('Enthalpy')
    title('He')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print('plot_N2','-dpng','-r0')
    
    % N2 PLOT
    figure
    plot(1:HX_slices, T(:, (HX_slices + 1 : 2 * HX_slices)))
    xlabel('Slices')
    ylabel('Enthalpy')
    title('N_2')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print('plot_HE','-dpng','-r0')
    
    % T_w PLOT
    figure
    plot(1:HX_slices, T_sol)
    xlabel('Slices')
    ylabel('Temperature')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print('plot_Twall','-dpng','-r0')

    
    % time iteration
	function [x, fval, exitflag] = heateq(j)  
		options = optimset('TolX', 1e-7, 'TolFun', 1e-7, ...
		    			   'MaxFunEvals', 1e7, 'MaxIter', 1e7, ...
		    			   'Display', 'iter');
        h_a_prev = h_sol(j - 1, 1 : HX_slices);
        h_b_prev = h_sol(j - 1, HX_slices + 1 : 2 * HX_slices);
        T_w_prev = T_sol(j - 1, :); 
        sol_guess = [h_a_prev h_b_prev T_w_prev]; 
		[x, fval, exitflag] = fsolve(@eqgen, sol_guess, options);

        % compute difference between the equation and zero
	    function F = eqgen(sol)
            F_a = zeros(1, HX_slices); 
            F_b = zeros(1, HX_slices); 
            F_w = zeros(1, HX_slices); 
            u_a = zeros(1, HX_slices); 
            u_b = zeros(1, HX_slices);
            u_a_prev = zeros(1, HX_slices); 
            u_b_prev = zeros(1, HX_slices);            
            h_a_delta = zeros(1, HX_slices); 
            h_b_delta = zeros(1, HX_slices); 
            
			h_a = sol(1 : HX_slices);
			h_b = sol(HX_slices + 1 : 2 * HX_slices);
            T_w = sol(2 * HX_slices + 1 : 3 * HX_slices);
            p_a = p_a_0;
            p_b = p_b_0;

			% energy equations
			for i = 1 : HX_slices  % slicing loop

				T_a = rp_thp(h_a(i), p_a(i), fluid_a);
				T_b = rp_thp(h_b(i), p_b(i), fluid_b);
             
				T_aw_delta = T_a - T_w;
				T_bw_delta = T_b - T_w;
				
				Q_cond_aw = HX_UA_aw / HX_slices * T_aw_delta;
				u_a(i) = rp_uhp(h_a(i), p_a(i), fluid_a);
				u_a_prev(i) = rp_uhp(h_a_prev(i), p_a(i), fluid_a);

				Q_cond_bw = HX_UA_bw / HX_slices * T_bw_delta;
				u_b(i) = rp_uhp(h_b(i), p_b(i), fluid_b);
				u_b_prev(i) = rp_uhp(h_b_prev(i), p_b(i), fluid_b);
                
                cp_w = cp(T_w); 
 
    			if i == 1
					h_a_delta(i) = h_a_in - h_a(1);
					h_b_delta(HX_slices) = h_b_in - h_b(HX_slices);
	            else
    				h_a_delta(i) = h_a(i - 1) - h_a(i);
    				h_b_delta(HX_slices + 1 - i) = h_b(HX_slices + 2 - i) - h_b(HX_slices + 1 - i);
	            end
            end
            
			for i = 1 : HX_slices
				F_a(i) = -(u_a(i) - u_a_prev(i)) * M / t_delta + (m * h_a_delta(i) + Q_cond_aw(i));
				F_b(i) = -(u_b(i) - u_b_prev(i)) * M / t_delta + (m * h_b_delta(i) + Q_cond_bw(i));
                F_w(i) = -(T_w(i) - T_w_prev(i)) * M_w / t_delta + (Q_cond_aw(i) + Q_cond_bw(i)) / cp_w(i) ;
            end          
                   
			% final exit vector
			F = [F_a  F_b  F_w];
		end
	end
end  % RN_01
