function [Q, data, T_a_sol, T_b_sol, T_w_sol] = RN_08_a2 
    % With radiative heat transfer 
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
    SWITCH = 0; % 1 for co-current, 0 for counter-current 
    t_delta = .1;  % time step, s
    t = 20;  % number of time steps, -
    HX_slices = 20;  % number of slices, -
    Wall_slices = 20;  % number of wall slices, -
    N = 1 + Wall_slices + 1;  % total length of slices across HX, - 
    lib = 'CP';  % property library to use

    % HX DATA
    m = 1;  % mass flow, kg/s
    M = 1;  % mass of streams content within a cell, kg
    M_w = 1;  % mass of wall section, kg 
    b_x = 1;  % length of wall section, m 
    HX_UA_data = {'nitrogen', 2500; ...
                    'helium', 2450}; % HX coefficient, W/K
    
    % RADIATION HEAT TRANSFER
    sigma = 5.676e-8; % Stefan-Boltzmann constant, W/m^2 * K^4
    As = 1; % surface area that is "seen" in radiation heat transfer, m^2
    F = 1; % view factor 
                
    % INITIAL DATA
    fluid_a = 'helium';  % stream A fluid name
    p_a_in = 101325;  % inlet pressure of stream A, Pa
    q_a_in = .765;  % inlet temperature of stream A, K
    h_a_in = prop_hqp(q_a_in, p_a_in, fluid_a, lib);  % inlet enthalpy of stream A, J/kg
    fluid_b = 'nitrogen';  % stream B fluid name
    p_b_in = 101325;  % inlet pressure of stream B, Pa
    T_b_in = 100;  % inlet temperature of stream B, K
    h_b_in = prop_htp(T_b_in, p_b_in, fluid_b, lib);  % inlet enthalpy of stream B, J/kg
    T_w_init = 60;  % initial wall temperature, K
    T_ext_init = 300; % exterior temperature, K
    T_a_in = 4; % Used only as a lower T limit for y-axis on plots

    % INITIAL CONDITIONS TIMES LENGTH OF HX
    h_a_0 = h_a_in * ones(HX_slices, 1);
    h_b_0 = h_b_in * ones(HX_slices, 1);
    p_a_0 = p_a_in * ones(HX_slices, 1);
    p_b_0 = p_b_in * ones(HX_slices, 1);
    T_ext = T_ext_init * ones(HX_slices, 1); 
    T_w_0 = T_w_init * ones(HX_slices, Wall_slices);
	    
    % COMBINE INITIAL CONDITIONS
    data(:, :, 1) = mcomb(h_a_0, T_w_0, h_b_0);
    % The data matrix keeps the whole solution
    % i are the HX_slices
    % j are the Wall_slices
    % k are the time steps 
        
    for k = 2 : t + 1 
        % h SOLVER
        [sol, fval, exitflag] = solve_h(k);
		disp(['Time step completed ' num2str(k-1) ...
            ' Time ' num2str(toc/60) ' min '...
            ' Fval ' num2str(sum(fval)) ...
            ' Exit Flag ' num2str(exitflag)])
        data(:, 1, k) = sol(:, 1); % stack 2-d to 3-d (stream A)
        data(:, N, k) = sol(:, 2); % stack 2-d to 3-d (stream B)
        
        % T SOLVER
        for i = 1 : HX_slices
            sol = solve_T(i,k);
            data(i, 2 : N - 1, k) = sol(end, :); % stack 
        end   
        
        % dump and re-load coolprop
        if rem(k,CP_dump) == 0
        unloadlibrary 'coolprop'
        loadcoolprop; 
        end 
    end
    
    % CONVERT TO T, Q & PLOT
    [T_a_sol, T_b_sol, T_w_sol, Q] = TQ_calc(data);   
    plots; 
    
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
    
    % FUNCTION THAT PULLS HX_UA DATA
    function value = HX_UA(fluid)        
        f = strcmp(fluid, HX_UA_data);   
        value = HX_UA_data{f, 2};
    end

    % FUNCTION THAT CALCULATES Q_cond 
    function Q_cond = Q_cond(h, p, T_w, fluid)
        T = propsc_thp(h, p, fluid, lib);
        T_delta = T_w - T;
        Q_cond = HX_UA(fluid) / HX_slices * T_delta;
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
        Q_rad = As * F * sigma * (T_ext.^4 - T.^4);
    end

    % FUNCTION THAT CALCULATES dudt
    function dudt = dudt(h, h_prev, p, fluid)
		dudt = propsc_uhp(h, p, fluid, lib) - propsc_uhp(h_prev, p, fluid, lib);
    end

    % FUNCTION THAT CONVERTS h,T to Q
    function [T_a_sol, T_b_sol, T_w_sol, Q] = TQ_calc(data)
        Q = zeros(HX_slices, 5); 
        T_a_sol = zeros(HX_slices, t + 1);
        T_b_sol = zeros(HX_slices, t + 1);
        
        for k = 1 : t + 1
            Q(:, 1, k) = Q_rad1(data(:, 1, k), data(:, N, k), p_a_0, p_b_0); % Q_rad_ab
            Q(:, 2, k) = Q_rad2(data(:, 1, k), p_a_0, fluid_a); % Q_rad_ae
            Q(:, 3, k) = Q_rad2(data(:, N, k), p_b_0, fluid_b); % Q_rad_be
            Q(:, 4, k) = Q_cond(data(:, 1, k), p_a_0, data(:, 2, k), fluid_a); % Q_cond_aw
            Q(:, 5, k) = Q_cond(data(:, N, k), p_b_0, data(:, N-1, k), fluid_b); % Q_cond_bw
            T_a_sol(:, k) = props_thp(data(:, 1, k), p_a_0, fluid_a, lib);  % j = 1
            T_b_sol(:, k) = props_thp(data(:, N, k), p_b_0, fluid_b, lib);  % j = N
        end
        
        T_w_sol = data(:, 2 : N - 1, :); % for j NOT than 1 or N 
    end 
 
    % ************************ PART III PLOTS ***************************
    function plots 
    % HE PLOT 
    plot(1 : HX_slices, T_a_sol, 'r')
    xlabel('Slices')
    ylabel('Temperature')
    title('He')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_N2' num2str(HX_slices)],'-dpng','-r0')
    
    % N2 PLOT
    figure
    plot(1 : HX_slices, T_b_sol, 'b')
    xlabel('Slices')
    ylabel('Temperature')
    title('N_2')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_HE' num2str(HX_slices)],'-dpng','-r0')
    
    % OVERALL PLOT
    figure 
    hold on 
    plot(1 : HX_slices, T_a_sol, 'k')
    plot(1 : HX_slices, T_b_sol, 'r')
    plot(1 : HX_slices, squeeze(T_w_sol(:, 1, :)), 'b')
    plot(1 : HX_slices, squeeze(T_w_sol(:, Wall_slices/2, :)), 'g')
    plot(1 : HX_slices, squeeze(T_w_sol(:, Wall_slices, :)), 'm')
    xlabel('Slices')
    ylabel('Temperature')
    title('Overall')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_overall' num2str(HX_slices)],'-dpng','-r0')
    
    % Q OVERALL PLOT
    figure
    hold on
    plot(1 : HX_slices, squeeze(Q(:, 1, :)), 'g') % Q_rad_ab
    plot(1 : HX_slices, squeeze(Q(:, 2, :)), 'm') % Q_rad_ae
    plot(1 : HX_slices, squeeze(Q(:, 3, :)), 'k') % Q_rad_be
    plot(1 : HX_slices, squeeze(Q(:, 4, :)), 'r') % Q_cond_aw  
    plot(1 : HX_slices, squeeze(Q(:, 5, :)), 'b') % Q_cond_bw
    xlabel('Slices')
    ylabel('Heat')
    title('Q Plot')
     
    % Q ZOOMED IN PLOT
    figure
    hold on
    plot(1 : HX_slices, squeeze(Q(:, 2, :)), 'm') % Q_rad_ae
    plot(1 : HX_slices, squeeze(Q(:, 3, :)), 'k') % Q_rad_be
    xlabel('Slices')
    ylabel('Heat')
    title('Q Plot')
    
    figure
    plot(1 : Wall_slices, squeeze(T_w_sol(HX_slices/2, :, :)), 'b')
    xlabel('Wall Slices')
    ylabel('Temperature')
    title('Wall Middle Plot')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(['plot_wall' num2str(HX_slices)],'-dpng','-r0')
    
    % 3D PLOT
    figure 
    [X,Y] = meshgrid(1 : HX_slices, 1 : N);
    time = [2, round(t/2), round(3 * t/4), t + 1]; % which time to plot
    colormap autumn 
    
    for I = 1:4 %4 sub-plots 
        subplot(2,2,I)
        surf(X,Y,[T_a_sol(:,time(I)) T_w_sol(:,:,time(I)) T_b_sol(:,time(I))]')
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
    
    % ************************ PART IV SOLVERS **************************
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
	    function eps = eqgen(sol)
            % pre-allocate 
			dhdx_a = zeros(HX_slices, 1);
			dhdx_b = zeros(HX_slices, 1);
            
            % split solution vector 
			h_a = sol(:, 1);
			h_b = sol(:, 2);
            p_a = p_a_0;
            p_b = p_b_0;
            T_w_a = T_w_prev(:, 1);
            T_w_b = T_w_prev(:, Wall_slices);

            dudt_a = dudt(h_a, h_a_prev, p_a, fluid_a);
            dudt_b = dudt(h_b, h_b_prev, p_b, fluid_b);
            Q_cond_aw = Q_cond(h_a, p_a, T_w_a, fluid_a);
            Q_cond_bw = Q_cond(h_b, p_b, T_w_b, fluid_b);
            Q_rad_ab = Q_rad1(h_a, h_b, p_a, p_b);
            Q_rad_ba = - Q_rad_ab; 
            Q_rad_ae = Q_rad2(h_a, p_a, fluid_a);
            Q_rad_be = Q_rad2(h_b, p_b, fluid_b);
            
			% enthalpy deltas - stream A
			dhdx_a(1) = h_a_in - h_a(1);
    		dhdx_a(2 : HX_slices) = h_a(1 : HX_slices - 1) - h_a(2 : HX_slices);
            
            % enthalpy deltas - stream B 
            if SWITCH == 1
                dhdx_b(1) = h_b_in - h_b(1);
                dhdx_b(2 : HX_slices) = h_b(1 : HX_slices - 1) - h_b(2 : HX_slices);
            elseif SWITCH == 0 
                dhdx_b(HX_slices) = h_b_in - h_b(HX_slices);
                dhdx_b(HX_slices - 1 : -1 : 1) = h_b(HX_slices : -1 : 2)...
                    - h_b(HX_slices - 1 : -1 : 1);
            end 
            
            % final equations 
			eps_a = - dudt_a * M / t_delta + m * dhdx_a ...
                    + Q_cond_aw + Q_rad_ab + Q_rad_ae;
			eps_b = - dudt_b * M / t_delta + m * dhdx_b ...
                    + Q_cond_bw + Q_rad_ba + Q_rad_be;
            
			% combine final exit vector
			eps = [eps_a  eps_b];
		end
    end

	function sol = solve_T(i,k)  % T solver                         
        
        % get data 
        W = Wall_slices; % short hand notation 
        [~, T_w_prev, ~] = msplit(data, k-1); 
        [h_a, ~, h_b] = msplit(data, k);
        p_a = p_a_0;
        p_b = p_b_0;
        
        % calculate Q from h 
        Q_cond_aw = Q_cond(h_a(i), p_a(i), T_w_prev(i, 1), fluid_a);   
        Q_cond_bw = Q_cond(h_b(i), p_b(i), T_w_prev(i, W), fluid_b);       
 
        % launch ode45 
        sol_guess = T_w_prev(i,:)';
        t_span = [0, t_delta];
        [~,sol] = ode45(@eqgen,t_span,sol_guess);

        % generate equations 
	    function dTdt_w = eqgen(~,T_w)
            
            % ensure we have column vectors 
            dTdx_w = zeros(Wall_slices, 1);
            Q_cond = zeros(Wall_slices, 1); 
            
            % Calculate K_w and cp
            K_w = K(T_w) * 100;
            cp_w = cp(T_w);
                 
            % Q_cond AT WALL EDGES 
            Q_cond(1) =  K_w(1)/b_x * (T_w(2) - T_w(1)) - Q_cond_aw;
            Q_cond(W) =  K_w(W)/b_x * (T_w(W-1) - T_w(W)) - Q_cond_bw;
                                    
            % Q_cond IN THE WALL
            dTdx_w(2 : W - 1) = T_w(3 : W) + T_w(1 : W - 2) - 2 * T_w(2 : W - 1);
            Q_cond(2 : W - 1) = K_w(2 : W - 1)/b_x .* dTdx_w(2 : W - 1);
            
            % Generate ODEs
            Q_rad = As * F * sigma * (T_ext.^4 - T_w.^4);
            dTdt_w = (Q_cond + Q_rad) ./ (M_w * cp_w);
        end
    end   
end  % RN_08_a2
