function [T_a_sol, T_b_sol, T_w_sol] = RN_06_a %With temperature distribution in the wall 
    addpath('H:\coolprop_lib')
    clc; clear;
    close all;
    tic
    
    % ************************ PART I DATA ********************************
    % INTEGRATOR DATA
    t_delta = 0.1;  % time step, s
    t = 20;  % number of time steps, -
    HX_slices = 10;  % number of slices, -
    Wall_slices = 1; % number of wall slices, -
    N = 1 + Wall_slices + 1; % total length of slices across HX, - 
    lib = 'CP'; % property library to use

    % HX DATA
    m = 1;  % mass flow, kg/s
    M = 1;  % mass of streams content within a cell, kg
    M_w = 1; % mass of wall section, kg 
    b_x = 1; % length of wall section, m 
    HX_UA_a = 5000; % HX coefficient, W/K
    HX_UA_b = 5000; % HX coefficient, W/K
    k = 5000; % wall conductance coefficient, W/K

    % INITIAL DATA
    fluid_a = 'helium';  % stream A fluid name
    p_a_in = 101325;  % inlet pressure of stream A, Pa
    T_a_in = 100;  % inlet temperature of stream A, K
    h_a_in = prop_htp(T_a_in, p_a_in, fluid_a, lib);  % inlet enthalpy of stream A, J/kg
    fluid_b = 'nitrogen';  % stream B fluid name
    p_b_in = 101325;  % inlet pressure of stream B, Pa
    T_b_in = 200;  % inlet temperature of stream B, K
    h_b_in = prop_htp(T_b_in, p_b_in, fluid_b, lib);  % inlet enthalpy of stream B, J/kg
    T_w_in = 150; % inlet wall temperature, K 

    % INITIAL CONDITIONS TIMES LENGTH OF HX 
    h_a_0 = h_a_in * ones(HX_slices, 1);
    h_b_0 = h_b_in * ones(HX_slices, 1);
    p_a_0 = p_a_in * ones(HX_slices, 1);
    p_b_0 = p_b_in * ones(HX_slices, 1);   
    T_w_0 = T_w_in * ones(HX_slices, Wall_slices); 
	    
    % COMBINE INITIAL CONDITIONS
    data(:, :, 1) = mcomb(h_a_0, T_w_0, h_b_0);
    % The data matrix keeps the whole solution
    % i are the HX_slices
    % j are the Wall_slices
    % k are the time steps 

    % *********************************************************************
    % CALL SOLVER
	for k = 2 : t + 1
		[sol, fval, exitflag] = heateq(k);
		disp(['Time step completed ' num2str(k-1)...
            ' Time ' num2str(toc/60) ' min '...
            ' Fval ' num2str(sum(fval)) ...
            ' Exit Flag ' num2str(exitflag)])
        data(:, :, k) = sol;
    end
    
    % GET FINAL SOLUTION AS T
    [T_a_sol, T_b_sol, T_w_sol] = Tconvert(data);
    
    
    % ************************ PART II FUNCTIONS **************************
    % FUNCTION THAT SPLITS MATRIX
    function [h_a, T_w, h_b] = msplit(matrix, k)
        h_a = matrix(:, 1, k);
        T_w = matrix(:, 2 : N - 1, k);
        h_b = matrix(:, N, k); 
    end

    % FUNCTION THAT COMBINES MATRICES
    function matrix = mcomb(A, B, C)
        matrix = [A B C];
    end    
    
    % FUNCTION THAT CALCULATES Q
    function [Q_cond, cp_w, T_w_delta] = ...
            QCpTw_calc(h_a, h_b, p_a, p_b, T_w, T_w_prev)   
        
        % pre-allocate 
        T = zeros(HX_slices, N);
        T_delta = zeros(HX_slices, N); 
        T_prev = zeros(HX_slices, N);
        T_w_delta = zeros(HX_slices, N);
        Q_cond = zeros(HX_slices, N);
        cp_w = zeros(HX_slices, N); 
                
        for i = 1 : HX_slices
            
            % Create a temperature matrix 
            T(i,1) = prop_thp(h_a(i), p_a(i), fluid_a, lib);
            T(i,N) = prop_thp(h_b(i), p_b(i), fluid_b, lib);
            T(i,2:N-1) = T_w(i,:);
            
            % Grab temperatures value from previous time step
            % (but only in the wall)
            T_prev(i,2:N-1) = T_w_prev(i,:); 
            
            % Q_cond AT EDGES
            T_delta(i,1) = T(i,2)-T(i,1);
            T_delta(i,N) = T(i,N-1)-T(i,N);
            Q_cond(i,1) = HX_UA_a / HX_slices * T_delta(i,1);
            Q_cond(i,N) = HX_UA_b / HX_slices * T_delta(i,N);          
            
            % Q_cond IN THE WALL 
            for j = 2 : N - 1
                T_delta(i,j) = 2*T(i,j) - T(i,j+1) - T(i,j-1);
                Q_cond(i,j) = k/b_x * T_delta(i,j);
            end
            
            % Calculate T_w_delta (in time) and evaluate C_p            
            for j = 2 : N - 1             
                T_w_delta(i,j) = T(i,j) - T_prev(i,j); 
                cp_w(i,j) = cp(T(i,j));
            end           
        end      
    end 

    % FUNCTION THAT COMPUTES delta_u
    function delta_u = u_delta(h, h_prev, p, fluid)
		delta_u = prop_uhp(h, p, fluid, lib) - prop_uhp(h_prev, p, fluid, lib); 
    end

    % FUNCTION THAT CONVERTS h TO T
    function [T_a_sol, T_b_sol, T_w_sol] = Tconvert(data)
        % pre-allocate 
        T_a_sol = zeros(HX_slices, t + 1);
        T_b_sol = zeros(HX_slices, t + 1);
    
        % convert
        T_w_sol = data(:, 2 : N - 1, :);
        for k = 1 : t + 1
            for i =  1 : HX_slices
                T_a_sol(i,k) = prop_thp(data(i,1,k), p_a_in, fluid_a, lib);
                T_b_sol(i,k) = prop_thp(data(i,N,k), p_b_in, fluid_b, lib);
            end
        end
    end
 
    % ************************ PART III PLOTS **************************
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
    plot(1:HX_slices, T_a_sol, 'k')
    plot(1:HX_slices, T_b_sol, 'b')
    plot(1:HX_slices, squeeze(T_w_sol(:, 1, :)), 'g')
    %plot(1:HX_slices, squeeze(T_w_sol(:, 2, :)), 'r')
    xlabel('Slices')
    ylabel('Temperature')
    title('Overall')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_overall' num2str(HX_slices)],'-dpng','-r0')
    
    
    % ************************ PART IV SOLVER **************************
	function [x, fval, exitflag] = heateq(k)  % time iteration
		options = optimset('TolX', 1e-7, 'TolFun', 1e-7, ...
		    			   'MaxFunEvals', 3000, 'MaxIter', 3000, ...
		    			   'Display', 'iter');
                       
        [h_a_prev, T_w_prev, h_b_prev] = msplit(data, k-1);
        sol_guess = mcomb(h_a_prev, T_w_prev, h_b_prev);
		[x, fval, exitflag] = fsolve(@eqgen, sol_guess, options);

        % compute difference between the equation and zero
	    function F = eqgen(sol)
            
            % pre-allocate 
   	    	F = zeros(HX_slices, N);
	    	u_a_delta = zeros(HX_slices, 1);
	    	u_b_delta = zeros(HX_slices, 1);
			h_a_delta = zeros(HX_slices, 1);
			h_b_delta = zeros(HX_slices, 1);

            % split solution vector 
			[h_a, T_w, h_b] = msplit(sol,1);
            p_a = p_a_0;
            p_b = p_b_0;
            
            % calculate Q_cond, cp_w, T_w_delta
            [Q_cond, cp_w, T_w_delta] = ...
                QCpTw_calc(h_a, h_b, p_a, p_b, T_w, T_w_prev);
                        
			% internal energy deltas 
			for i = 1 : HX_slices  % slicing loop                
                u_a_delta(i) = u_delta(h_a(i), h_a_prev(i), p_a(i), fluid_a);
                u_b_delta(i) = u_delta(h_b(i), h_b_prev(i), p_b(i), fluid_b);
                              
            % enthalpy deltas  
    			if i == 1
					h_a_delta(1) = h_a_in - h_a(1);
					h_b_delta(HX_slices) = h_b_in - h_b(HX_slices);
	            else
    				h_a_delta(i) = h_a(i - 1) - h_a(i);
    				h_b_delta(HX_slices + 1 - i) = ...
                        h_b(HX_slices + 2 - i) - h_b(HX_slices + 1 - i);
	            end
            end
            
            % final equations 
			for i = 1 : HX_slices % slicing loop  
                
				F(i,1) = - u_a_delta(i) * M / t_delta + m * h_a_delta(i) + Q_cond(i,1);
				F(i,N) = - u_b_delta(i) * M / t_delta + m * h_b_delta(i) + Q_cond(i,N);
                
                for j = 2 : N - 1 % wall temperatures in the middle 
                    F(i,j) = - T_w_delta(i,j) * M_w / t_delta + Q_cond(i,j) / cp_w(i,j) ;
                end
            end   
            
		end
	end
end  % RN_06_a
