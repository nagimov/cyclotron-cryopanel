function [data, p_a_data, p_b_data, T_a_sol, T_b_sol, T_w_sol, T_w3_data, Q] = RN_09_b 
    % Wall tail - with moving node
    % For low T, use HX_UA = 1500 and a factor of 100 for K, p_atm/3
    clc; clear;
    close all;
    tic
    
    % ************************** COOL PROP ********************************
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
    
    % ************************ PART I DATA ********************************
    % INTEGRATOR DATA
    t_delta = .1;  % time step, s
    t = 20;  % number of time steps, -
    HX_slices = 20;  % number of slices, -
    Wall_slices_1 = 10;  % number of wall slices, A to node, -
    Wall_slices_2 = 10;  % number of wall slices, B to node, - 
    Wall_slices_3 = 10;  % number of wall slices, node to tail, -
    W1 = Wall_slices_1; % running sum & short hand notation 
    W2 = W1 + Wall_slices_2; % running sum & short hand notation 
    W3 = W2 + Wall_slices_3; % running sum & short hand notation 
    N = 1 + W2 + 1;  % total number of slices across HX, -
    lib = 'CP';  % property library to use

    % HX DATA
    m = 1;  % mass flow, kg/s
    M = 1;  % mass of streams content within a cell, kg
    M_w = 1;  % mass of wall section, kg 
    b_x = 1;  % length of wall section, m 
    
    % RADIATION HEAT TRANSFER
    sigma = 5.676e-8; % Stefan-Boltzmann constant, W/m^2 * K^4
    
    % RADIATION HEAT TRANSFER - WALL
    rad_data_w = [1, 1;...
                  1, 1;...
                  10, 1];
    % column 1 is As
    % column 2 is F
    % rows are wall sections 
    As_w = built_rad_data(rad_data_w, 1);
    F_w = built_rad_data(rad_data_w, 2);
    
    % RADIATION HEAT TRANSFER - STREAMS 
    As = 1; % surface area that is "seen" in radiation heat transfer, m^2
    F = 1; % view factor 
                
    % INITIAL DATA
    p_atm = 101325; % define atmospheric pressure to be used as a reference, Pa 
    fluid_a = 'helium';  % stream A fluid name
    p_a_in = 1 * p_atm;  % inlet pressure of stream A, Pa
    T_a_in = 100;  % inlet temperature of stream A, K
    h_a_in = prop_htp(T_a_in, p_a_in, fluid_a, lib);  % inlet enthalpy of stream A, J/kg
    fluid_b = 'nitrogen';  % stream B fluid name
    p_b_in = 1 * p_atm;  % inlet pressure of stream B, Pa
    T_b_in = 200;  % inlet temperature of stream B, K
    h_b_in = prop_htp(T_b_in, p_b_in, fluid_b, lib);  % inlet enthalpy of stream B, J/kg
    T_w_init = 150;  % initial wall temperature, K
    T_ext_init = 300; % exterior temperature, K
    
    % SET NOMINAL VALUES
    m_nom = 1; % nominal mass flow rate, kg/s
    T_nom = 150; % nomial temperature, K            *** need to use 'q' 
    p_nom = p_atm; % nomial pressure, Pa
    HX_UA_nom_data = {'nitrogen', 3000; ...
                        'helium', 3000}; % nomial HX coefficient, W/K
    delta_p_nom_data = {'helium', p_atm / 3; ... 
                        'nitrogen', p_atm / 3}; % nomial pressure drop, Pa
    nom_values = {'helium', nom_calc(fluid_a);...
                   'nitrogen', nom_calc(fluid_b)}; % calculate other nominal values
    
    % INITIAL CONDITIONS TIMES LENGTH OF HX
    h_a_0 = h_a_in * ones(HX_slices, 1);
    h_b_0 = h_b_in * ones(HX_slices, 1);
    p_a_0 = p_dist(delta_p_nom_data{1,2}, p_a_in);
    p_b_0 = p_dist(delta_p_nom_data{2,2}, p_b_in);
    T_ext = T_ext_init * ones(max(HX_slices, W3), 1); 
    T_w1_0 = T_w_init * ones(HX_slices, Wall_slices_1);
    T_w2_0 = T_w_init * ones(HX_slices, Wall_slices_2);
    T_w3_0 = T_w_init * ones(HX_slices, Wall_slices_3);
	    
    % COMBINE INITIAL CONDITIONS FOR h and T
    data(:, :, 1) = [h_a_0, T_w1_0, T_w2_0, h_b_0];
    T_w3_data(:, :, 1) = T_w3_0; 
    % The data matrix keeps the whole solution
    % i are the HX_slices
    % j are the Wall_slices
    % k are the time steps 
    
    % INTIAL CONDITIONS FOR p
    p_a_data(:, 1) = p_a_0;
    p_b_data(:, 1) = p_b_0; 
    % rows the HX_slices
    % last row is the delta_p
    % columns the time steps 
        
    for k = 2 : t + 1
       
        % p SOLVER
        sol = solve_p(k);
        p_a_data(:, k) = sol(:, 1); % stack 2-d to 3-d in time (stream A)
        p_b_data(:, k) = sol(:, 2); % stack 2-d to 3-d in time (stream B)
        
        % h SOLVER
        [sol, fval, exitflag] = solve_h(k);
		disp(['Time step completed ' num2str(k-1) ...
            ' Time ' num2str(toc/60) ' min '...
            ' Fval ' num2str(sum(fval)) ...
            ' Exit Flag ' num2str(exitflag)])
        data(:, 1, k) = sol(1 : 1 * HX_slices); % stack 2-d to 3-d (stream A)
        data(:, N, k) = sol(1 * HX_slices + 1 : 2 * HX_slices); % stack 2-d to 3-d (stream B)
        
        % T SOLVER
        for i = 1 : HX_slices
            [time,sol] = solve_T(i,k);
            data(i, 2 : N - 1, k) = sol(end, 1 : W2); % stack 2-d to 3-d in time (main wall)
            T_w3_data(i, :, k) = sol(end, W2 + 1 : W3); % stack 2-d to 3-d in time (wall tail)
        end
        
        % dump and re-load coolprop
        if rem(k,CP_dump) == 0 % check to see if time step is a multiple of CP_dump
        unloadlibrary 'coolprop'
        loadcoolprop; 
        end 
    end
    
    % CONVERT TO T, Q & PLOT
    [T_a_sol, T_b_sol, T_w_sol, Q] = TQ_calc(data);   
    plots_T; 
    plots_p; 
    
    % ************************ PART II FUNCTIONS **************************
    % FUNCTION THAT SPLITS MATRIX
    function [h_a, T_w, h_b] = msplit(matrix, k)
        h_a = matrix(:, 1, k);
        T_w = matrix(:, 2 : N - 1, k);
        h_b = matrix(:, N, k); 
    end

    % FUNCTION THAT CONVERTS RAD DATA INTO PROPER FORMAT (for the wall)
    function value = built_rad_data(rad_data_w, no)
        value = [rad_data_w(1, no) * ones(Wall_slices_1, 1);...
                rad_data_w(2, no) * ones(Wall_slices_2, 1);...
                rad_data_w(3, no) * ones(Wall_slices_3, 1)];
    end
    
    % FUNCTION THAT CALCULATES NOMINAL DATA 
    function values = nom_calc(fluid)        
        f = strcmp(fluid, HX_UA_nom_data);
        HX_UA_nom = HX_UA_nom_data{f, 2};
        delta_p_nom = delta_p_nom_data{f, 2};
        h_nom = prop_htp(T_nom, p_nom, fluid, lib); 
        [L_nom, mu_nom, rho_nom, Pr_nom] = ...
            propsc_LmurhoPr_hp(h_nom, p_nom, fluid, lib);
        values = [HX_UA_nom, delta_p_nom, L_nom, mu_nom, rho_nom, Pr_nom]; 
    end

    % FUCNTION THAT CALCULATES HX_UA
    function value = HX_UA(h, p, fluid)
        
        % grab nomimal values and calculate current values 
        f = strcmp(fluid, nom_values);           
        HX_UA_nom = nom_values{f, 2}(1);
        L_nom = nom_values{f, 2}(3);
        mu_nom = nom_values{f, 2}(4);
        Pr_nom = nom_values{f, 2}(6);
        [L,mu,~,Pr] = propsc_LmurhoPr_hp(h, p, fluid, lib);
        
        % scale HX_UA
        value = HX_UA_nom * (L / L_nom) .* (m / m_nom).^.8 .*...
            (mu_nom ./ mu).^.8 .* (Pr / Pr_nom).^(1/3);
    end

    % FUNCTION THAT CALCULATES Q_cond 
    function Q_cond = Q_cond(h, p, T_w, fluid)
        T = propsc_thp(h, p, fluid, lib);
        T_delta = T_w - T;
        Q_cond = HX_UA(h, p, fluid) / HX_slices .* T_delta;
    end 

    % FUNCTION THAT CALCULATES Q_rad (between two streams)
    function Q_rad = Q_rad1(h_a, h_b, p_a, p_b)
        T_a = propsc_thp(h_a, p_a, fluid_a, lib);
        T_b = propsc_thp(h_b, p_b, fluid_b, lib);
        Q_rad = As * F * sigma * (T_b.^4 - T_a.^4);
    end

    % FUNCTION THAT CALCULATES Q_rad (between streams and an exterior) 
    function Q_rad = Q_rad2(h, p, fluid)
        T = propsc_thp(h, p, fluid, lib);
        Q_rad = As * F * sigma * (T_ext(1: HX_slices).^4 - T.^4);
    end

    % FUNCTION THAT CALCULATES dudt
    function dudt = dudt(h, h_prev, p, fluid)
		dudt = propsc_uhp(h, p, fluid, lib) - propsc_uhp(h_prev, p, fluid, lib);
    end

    % FUNCTION THAT CALCULATES delta_p
    function delta_p = p_calc(h, p, p_in, fluid)
               
        % grab nomimal vlaues and calculate current values 
        f = strcmp(fluid, nom_values);
        delta_p_nom = nom_values{f, 2}(2);
        mu_nom = nom_values{f, 2}(4);
        rho_nom = nom_values{f, 2}(5);
        h_mean = mean(h);
        p_mean = mean(p(1:HX_slices));
        [~,mu,rho,~] = propsc_LmurhoPr_hp(h_mean, p_in, fluid, lib);
               
        % calculate delta_p
        delta_p = delta_p_nom * (m / m_nom)^1.8 * (mu / mu_nom)^0.2 ...
            * (rho_nom / rho);
    end 

    % FUNCTION THAT CALCULATES p_dist 
    function p_dist = p_dist(delta_p, p_in)
        delta_p = delta_p / HX_slices; % delta_p per slice 
        p = p_in : -delta_p : p_in - delta_p * (HX_slices - 1);
        p_dist = [p, delta_p * HX_slices]; 
    end

    % FUNCTION THAT CONVERTS h,T to Q
    function [T_a_sol, T_b_sol, T_w_sol, Q] = TQ_calc(data)
        Q = zeros(HX_slices, 5); 
        T_a_sol = zeros(HX_slices, t + 1);
        T_b_sol = zeros(HX_slices, t + 1);
        T_w_sol = data(:, 2 : N - 1, :); % for j NOT than 1 or N 
        
        for k = 1 : t + 1
            T_a_sol(:,k) = propsc_thp(data(:, 1, k), p_a_0, fluid_a, lib);  % T for j = 1
            T_b_sol(:,k) = propsc_thp(data(:, N, k), p_b_0, fluid_b, lib);  % T for j = N
        end
        
    end 
 
    % ************************ PART III PLOTS ***************************
    function plots_T 
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
    
    % OVERALL PLOT
    figure 
    hold on 
    plot(1:HX_slices, T_a_sol, 'k')
    plot(1:HX_slices, T_b_sol, 'r')
    plot(1:HX_slices, squeeze(T_w_sol(:, 1, :)), 'b')
    plot(1:HX_slices, squeeze(T_w_sol(:, W1, :)), 'g')
    plot(1:HX_slices, squeeze(T_w_sol(:, W2, :)), 'm')
    xlabel('Slices')
    ylabel('Temperature')
    title('Overall')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_overall' num2str(HX_slices)],'-dpng','-r0')
    
    % WALL TAIL PLOT
    figure
    hold on
    plot(W2 + 1 : W3, squeeze(T_w3_data(HX_slices/2, :, :)),'g')
    plot(W2 + 1 : W3, squeeze(T_w3_data(HX_slices/2, :, :)),'r')
    xlabel('Wall Slices')
    ylabel('Temperature')
    title('Wall Tail Plot')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_tail' num2str(HX_slices)],'-dpng','-r0')
    
    % WALL CROSS SECTION PLOT 
    figure
    plot(1 : W2, squeeze(T_w_sol(HX_slices/2, :, :)), 'b')
    xlabel('Wall Slices')
    ylabel('Temperature')
    title('Wall Middle Plot')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_wall' num2str(HX_slices)],'-dpng','-r0')
    
    % 3D PLOT
    figure 
    [X,Y] = meshgrid(1 : HX_slices, 1 : N); % get a mesh grid
    time = [1+1, round(t/2), round(3*t/4), t+1]; % which times to plot
    colormap autumn 
    
    for I = 1:4 %4 sub-plots 
        subplot(2, 2, I)
        surf(X,Y,[T_a_sol(:,time(I)) T_w_sol(:, :, time(I)) T_b_sol(:, time(I))]')
        xlabel('Length Slices')
        ylabel('Wall Slices')
        zlabel('Temperature')
        title(['Time step ' num2str(time(I))])
        axis([1 HX_slices 1 N T_a_in T_b_in])
        shading interp;
        caxis manual
        caxis([T_a_in T_b_in]); % fix color distribution to always be between Tmin, Tmax
        colorbar 
    end 
    end
    
    function plots_p
    % p vs time PLOT
    figure
    hold on 
    plot(2 : t + 1, p_atm * ones(1, t),'k')
    plot(2 : t + 1, p_a_data(1 : HX_slices - 1, 2 : t + 1)','r')   
    plot(2 : t + 1, p_b_data(1 : HX_slices - 1, 2 : t + 1)','b')
    plot(2 : t + 1, p_a_data(HX_slices, 2 : t + 1)','r:')   
    plot(2 : t + 1, p_b_data(HX_slices, 2 : t + 1)','b:')
    xlabel('Time')
    ylabel('Pressure')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_p' num2str(HX_slices)],'-dpng','-r0')
    
    % delta_p PLOT
    figure
    hold on
    plot(2 : t + 1, p_a_data(HX_slices + 1, 2 : t + 1),'r')
    plot(2 : t + 1, p_b_data(HX_slices + 1, 2 : t + 1),'b')
    xlabel('Time')
    ylabel('Pressure Drop')
    title('delta p Plot')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_delta_p' num2str(HX_slices)],'-dpng','-r0')
    end
    
    % ************************ PART IV SOLVERS **************************
    function sol = solve_p(k) % p solver
        [h_a_prev, ~, h_b_prev] = msplit(data, k-1);
        p_a_prev = p_a_data(:, k-1);
        p_b_prev = p_b_data(:, k-1);
        delta_p_a = p_calc(h_a_prev, p_a_prev, p_a_in, fluid_a);
        delta_p_b = p_calc(h_b_prev, p_b_prev, p_b_in, fluid_b);
        sol = [p_dist(delta_p_a, p_a_in)', p_dist(delta_p_b, p_b_in)'];
    end
    
    function [x, fval, exitflag] = solve_h(k)  % h solver

        % pull previous data and guess 
        [h_a_prev, T_w_prev, h_b_prev] = msplit(data, k-1);
        sol_guess = [h_a_prev, h_b_prev];
        
        % launch solver
        options = optimset('TolX', 1e-7, 'TolFun', 1e-7, ...
		    			   'MaxFunEvals', 1e7, 'MaxIter', 1e7, ...
		    			   'Display', 'iter');
		[x, fval, exitflag] = fsolve(@eqgen, sol_guess, options);

        % compute difference between the equation and zero
	    function F = eqgen(sol)
            
            % pre-allocate 
			dhdx_a = zeros(HX_slices, 1);
			dhdx_b = zeros(HX_slices, 1);
            
            % split solution vector 
			h_a = sol(:, 1);
			h_b = sol(:, 2);
            
            % get p and T_w data             
            T_w_a = T_w_prev(:, 1);
            T_w_b = T_w_prev(:, W2);
            p_a = p_a_data(1 : HX_slices, k);
            p_b = p_b_data(1 : HX_slices, k);
            
            % dudt
            dudt_a = dudt(h_a, h_a_prev, p_a, fluid_a);
            dudt_b = dudt(h_b, h_b_prev, p_b, fluid_b);
            
            % Q_cond
            Q_cond_aw = Q_cond(h_a, p_a, T_w_a, fluid_a);
            Q_cond_bw = Q_cond(h_b, p_b, T_w_b, fluid_b);
            
            % Q_rad
            Q_rad_ab = Q_rad1(h_a, h_b, p_a, p_b);
            Q_rad_ba = - Q_rad_ab;
            Q_rad_ae = Q_rad2(h_a, p_a, fluid_a);
            Q_rad_be = Q_rad2(h_b, p_b, fluid_b);
            
			% enthalpy deltas
			for i = 1 : HX_slices  % slicing loop                
    			if i == 1
					dhdx_a(1) = h_a_in - h_a(1);
					dhdx_b(HX_slices) = h_b_in - h_b(HX_slices);
	            else
    				dhdx_a(i) = h_a(i - 1) - h_a(i);
    				dhdx_b(HX_slices + 1 - i) = h_b(HX_slices + 2 - i) - h_b(HX_slices + 1 - i);
	            end
            end
                       
            % final equations 
			F_a = - dudt_a * M / t_delta + m * dhdx_a ...
                    + Q_cond_aw + Q_rad_ab + Q_rad_ae;
			F_b = - dudt_b * M / t_delta + m * dhdx_b ...
                    + Q_cond_bw + Q_rad_ba + Q_rad_be;
            
			% combine final exit vector
			F = [F_a  F_b];
		end
    end

	function [time,sol] = solve_T(i, k)  % T solver, for a fixed HX_slice                          
        
        % get data  
        [~, T_w_prev, ~] = msplit(data, k-1);
        T_w3_prev = T_w3_data(:, :, k-1);
        [h_a, ~, h_b] = msplit(data, k);
        p_a = p_a_data(1 : HX_slices, k);
        p_b = p_b_data(1 : HX_slices, k);
               
        % calculate Q from h 
        Q_cond_aw = Q_cond(h_a(i), p_a(i), T_w_prev(i, 1), fluid_a);   
        Q_cond_bw = Q_cond(h_b(i), p_b(i), T_w_prev(i, W2), fluid_b);       
 
        % launch ode45 
        sol_guess = [T_w_prev(i, :)'; T_w3_prev(i, :)'];
        size(sol_guess);
        t_span = [0, t_delta];
        [time,sol] = ode45(@eqgen, t_span, sol_guess);

        % generate equations 
	    function dTdt_w = eqgen(~,T_w)
            
            % ensure we have column vectors 
            dTdx_w = zeros(W3, 1);
            Q_cond = zeros(W3, 1); 
            
            % Calculate K_w and cp
            K_w = K(T_w);
            cp_w = cp(T_w);
                 
            % Q_cond AT WALL EDGES 
            Q_cond(1) =  K_w(1)/b_x * (T_w(2) - T_w(1)) - Q_cond_aw;
            Q_cond(W2) = K_w(W2)/b_x * (T_w(W2-1) - T_w(W2)) - Q_cond_bw;
            
            % Q_cond AT THE SPLIT
            Q_cond(W2 + 1) = K_w(W2 + 1)/b_x * (T_w(W1) + T_w(W2 + 2) - 2 * T_w(W2 + 1));
            Q_cond(W1) = K_w(W1)/b_x ...
                        * (T_w(W2 + 1) + T_w(W1 - 1) + T_w(W1 + 1) - 3 * T_w(W1));
                    
            % Q_cond AT THE END OF THE TAIL 
            Q_cond(W3) = K_w(W3)/b_x * (T_w(W3 - 1) - T_w(W3)); 
                                    
            % Q_cond IN THE MAIN WALL 
            dTdx_w(2 : W1 - 1)     = T_w(3 : W1) + T_w(1 : W1 - 2) ...
                                    - 2 * T_w(2 : W1 - 1);
            dTdx_w(W1 + 1 : W2 - 1) = T_w(W1 + 2 : W2) + T_w(W1 : W2 - 2) ...
                                    - 2 * T_w(W1 + 1 : W2 - 1);
            Q_cond(2 : W1 - 1)     = K_w(2 : W1 - 1) / b_x ...
                                    .* dTdx_w(2 : W1 - 1);
            Q_cond(W1 + 1 : W2 - 1) = K_w(W1 + 1 : W2 - 1) / b_x ...
                                    .* dTdx_w(W1 + 1 : W2 - 1); 
            
            % Q_cond IN THE TAIL
            dTdx_w(W2 + 2 : W3 - 1) = T_w(W2 + 3 : W3) + T_w(W2 + 1 : W3 - 2) ...
                                      - 2 * T_w(W2 + 2 : W3 - 1);
            Q_cond(W2 + 2 : W3 - 1) = K_w(W2 + 2 : W3 - 1)/b_x ...
                                      .* dTdx_w(W2 + 2 : W3 - 1); 
            
            % Q_rad & Generate ODEs
            Q_rad = As_w .* F_w * sigma .* (T_ext(1: W3).^4 - T_w.^4);
            dTdt_w = (Q_cond + Q_rad) ./ (M_w * cp_w);
        end
    end   
end  % RN_09_b
