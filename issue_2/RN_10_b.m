function [data, T_data, T_v_data, p_data] = RN_10_b
    % First version of the 4-stream model 
    clc; clear;
    close all;
    tic
    
    % ************************** COOL PROP ********************************
    CP_dump = 20; % number of time steps before we dump & re-load the library 
    path_to_lib = 'D:\CoolProp_wrapper_fast'; % path to coolprop shared library
    path_to_include= 'D:\CoolProp_wrapper_fast'; % path to coolprop's folder
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
                path_to_include,'alias','coolprop'); 
                % loading library with alias coolprop
        end
    end
    
    % ************************ SOLVER, HX DATA ****************************
    % SOLVER DATA
    SWITCH = 0; % 1 for co-current, 0 for counter-current 
    t_delta = .1;  % time step, s
    t = 20;  % number of time steps, -
    hx_slices = [20, 20];  % number of slices, -
    wall_slices = [10, 10; ... % stream A(C) to node, -
                   10, 10; ... % stream node to stream B(D), -
                   10, 10]; % node to tail, -   
                % left, right sides 
   
    % CALCULATE WALL SLICE INTERMEDIATES 
    w1 = wall_slices(1, :); % from stream A(C) to node
    w2 = w1 + wall_slices(2, :); % from node to stream B(D)
    w3 = w2 + wall_slices(3, :); % from node to tail 
    n = 1 * ones(1, 2) + w2 + 1 * ones(1, 2); 
    % (overall number of slices across HX)

    % HX DATA
    m = 1;  % mass flow, kg/s
    M = 1;  % mass of streams content within a cell, kg
    M_w = 1;  % mass of wall section, kg 
    b_x = 1;  % length of wall section, m 
                
    % ************************ INITIAL DATA *******************************
    % dimension j & ii are introduced 
    p_atm = 101325; % define atmospheric pressure to be used as a reference, Pa 
    fluid = {'helium', 'helium';  % stream A, stream C
            'nitrogen', 'nitrogen'};  % stream B, stream D
    p_in = [1 * p_atm, 1 * p_atm;   % stream A, stream C; pressure, Pa
            1 * p_atm, 1 * p_atm];  % stream B, stream D; pressure, Pa
    T_in = [100, 50;               % stream A, stream C; temperature, K
            200, 100];             % stream B, stream D; temperature, K
    h_in = hcalc(T_in, p_in, fluid); % convert input to h     
    T_w_init = [150, 75]; % initial wall temperature (left, right), K, -   
    T_v_init = [150, 75]; % initial wall tail temperature (left, right), K, - 
    
    % ************************ RADIATION DATA *****************************
    T_ext_init = 300; % External body temperature, K 
    sigma = 5.676e-8; % Stefan-Boltzmann constant, W/m^2 * K^4
    As = [1, 1]; % Surface Area, for He streams only, m^2
    F = [1, 1]; % View factor, for He streams only, - 
    As_data = [1, 1;...
               1, 1;...
               10, 10];
    % Surface area that is "seen" in radiative heat transfer, m^2 
    % rows are wall sections
    % columns are left, right
    F_data = [1, 1;...
              1, 1;...
              1, 1];
    % Radiation view factor
    % rows are wall sections
    % columns are left, right
    
    % ************************ NOMINAL DATA *******************************
    m_nom = 1; % nominal mass flow rate, kg/s
    T_nom = 100; % nomial temperature, K            
    p_nom = p_atm; % nomial pressure, Pa
    HX_UA_nom_data = {'nitrogen', 3000; ...
                        'helium', 3000}; % nomial HX coefficient, W/K
    delta_p_nom_data = {'helium', p_atm / 3; ... 
                        'nitrogen', p_atm / 3}; % nomial pressure drop, Pa
    nom_values = {'helium', nom_calc('helium'); ...
                   'nitrogen', nom_calc('nitrogen')}; % calculate other nominal values
               
    % ************************ PRE-ALLOCATE *******************************
    % The data matrices keep the whole solution
    % i are the HX_slices
    % j are the Wall_slices
    % k are the time steps 
    % ii is the side of the cryo-panel 
    data = zeros(max(hx_slices), max(n), t + 1, 2);
    
    % data matrix converted into T
    T_data = zeros(max(hx_slices), max(n), t + 1, 2);
    
    % stores wall tail T
    T_v_data = zeros(max(hx_slices), max(w3-w2), t + 1, 2);
    
    % stores pressure 
    p_data = zeros(max(hx_slices), 2, t + 1, 2);
    h_0 = zeros(max(hx_slices), 2, 2);
    p_0 = zeros(max(hx_slices), 2, 2);

    % ************************ MAIN LOOPS *********************************
    for ii = 1 : 2
        
        % Pull the right number of slices & rad data  
        HX_slices = hx_slices(ii);
        N = n(ii);  W1 = w1(ii);  W2 = w2(ii);  W3 = w3(ii);
        As_w = built_rad_data(As_data(:, ii));
        F_w = built_rad_data(F_data(:, ii));
        
        % EXTEND INITIAL CONDITIONS TO ALL SLICES 
        % (dimension i is introduced)
        % (dimension j is extended)
        h_0(:, 1, ii) = h_in(1, ii) * ones(HX_slices, 1); 
        h_0(:, 2, ii) = h_in(2, ii) * ones(HX_slices, 1);
        p_0(:, 1, ii) = p_in(1, ii) * ones(HX_slices, 1);
        p_0(:, 2, ii) = p_in(2, ii) * ones(HX_slices, 1);
        T_w_0 = T_w_init(ii) * ones(HX_slices, W2);
        T_v_0 = T_v_init(ii) * ones(HX_slices, W3 - W2);
        T_ext = T_ext_init * ones(max(HX_slices, W3), 1); 
                
        % COMBINE INITIAL CONDITIONS 
        % (dimension k is introduced)
        data(:, :, 1, ii) = [h_0(:, 1, ii), T_w_0, h_0(:, 2, ii)];
        p_data(:, :, 1, ii) = p_0(:, : , ii); 
        T_v_data(:, :, 1, ii) = T_v_0; 
        
        for k = 2 : t + 1
            
            % h SOLVER
            [sol, fval, exitflag] = solve_h(k, ii);
            disp(['Time step completed ' num2str(k - 1) '(' num2str(ii) ')'...
                ' Time ' num2str(toc/60) ' min '...
                ' Fval ' num2str(sum(fval)) ...
                ' Exit Flag ' num2str(exitflag)])
            
            % stack in time 
            data(:, 1, k, ii) = sol(:, 1); 
            data(:, N, k, ii) = sol(:, 2); 
            
            % T SOLVER
            for i = 1 : HX_slices
                sol = solve_T(i, k, ii);
                
                %stack in time
                data(i, 2 : N - 1, k, ii) = sol(end, 1 : W2); % main wall
                T_v_data(i, :, k, ii) = sol(end, W2 + 1 : W3); % wall tail
            end
            
            % convert h to T (all time steps except k = 1)
            T_data(:, 1, k, ii) = props_thp(data(:, 1, k, ii), ...
                p_data(:, 1, 1, ii), fluid{1, ii});  % T for j = 1
            T_data(:, N, k, ii) = props_thp(data(:, N, k, ii), ...
                p_data(:, 1, 1, ii), fluid{2, ii});  % T for j = N
            
            % dump and re-load coolprop
            if rem(k, CP_dump) == 0 
                % check to see if time step is a multiple of CP_dump
                unloadlibrary 'coolprop'
                loadcoolprop; 
            end 
        end
        
        % get T for k = 1
        T_data(:, 1, 1, ii) = T_in(1, ii) * ones(HX_slices, 1);
        T_data(:, N, 1, ii) = T_in(2, ii) * ones(HX_slices, 1);
    end 
    
    % copy wall temperatures to T_data matrix
    T_data(:, 2 : N - 1, :, :) = data(:, 2 : N - 1, :, :); % for j NOT 
    
    % ***************************** FUNCTIONS *****************************
    % Enthalpy conversion  
    function h = hcalc(T, p, fluid)
        h = zeros(2); 
        for j = 1 : 2
            for ii = 1 : 2
                h(j, ii) = props_htp(T(j, ii), p(j, ii), fluid{j, ii});
            end
        end 
    end 

    % Split data matrix 
    function [h_a, T_w, h_b] = msplit(matrix, k, ii)
        h_a = matrix(:, 1, k, ii);
        T_w = matrix(:, 2 : N - 1, k, ii);
        h_b = matrix(:, N, k, ii); 
    end

    % Calculate nominal values  
    function values = nom_calc(fluid)        
        f = strcmp(fluid, HX_UA_nom_data);
        HX_UA_nom = HX_UA_nom_data{f, 2};
        delta_p_nom = delta_p_nom_data{f, 2};
        h_nom = prop_htp(T_nom, p_nom, fluid, 'CP'); 
        [L_nom, mu_nom, rho_nom, Pr_nom] = ...
            propsc_LmurhoPr_hp(h_nom, p_nom, fluid);
        values = [HX_UA_nom, delta_p_nom, L_nom, mu_nom, rho_nom, Pr_nom]; 
    end

    % Calculate HX_UA
    function value = HX_UA(h, p, fluid)
        
        % grab nomimal values and calculate current values 
        f = strcmp(fluid, nom_values);           
        HX_UA_nom = nom_values{f, 2}(1);
        L_nom = nom_values{f, 2}(3);
        mu_nom = nom_values{f, 2}(4);
        Pr_nom = nom_values{f, 2}(6);
        [L,mu,~,Pr] = propsc_LmurhoPr_hp(h, p, fluid);
        
        % scale HX_UA
        value = HX_UA_nom * (L / L_nom) .* (m / m_nom).^.8 .*...
            (mu_nom ./ mu).^.8 .* (Pr / Pr_nom).^(1/3);
    end

    % FUNCTION THAT CONVERTS RAD DATA INTO PROPER FORMAT (for the wall)
    function value = built_rad_data(rad_data_w)
        value = [rad_data_w(1) * ones(W1, 1);...
                rad_data_w(2) * ones(W2 - W1, 1);...
                rad_data_w(3) * ones(W3 - W2, 1)];
    end

    % FUNCTION THAT CALCULATES Q_cond 
    function Q_cond = Q_cond(h, p, T_w, fluid)
        T = propsc_thp(h, p, fluid);
        T_delta = T_w - T;
        Q_cond = HX_UA(h, p, fluid) / HX_slices .* T_delta;
    end

    % FUNCTION THAT CALCULATES Q_rad (between streams and an exterior) 
    function Q_rad = Q_rad2(h, p, fluid, ii)
        T = propsc_thp(h, p, fluid);
        Q_rad = As(ii) * F(ii) * sigma * (T_ext(1: HX_slices).^4 - T.^4);
    end

    % FUNCTION THAT CALCULATES dudt
    function dudt = dudt(h, h_prev, p, fluid)
		dudt = propsc_uhp(h, p, fluid) - propsc_uhp(h_prev, p, fluid);
    end

    %  Pull title for plots 
    function value = pull_title(ii)
        if ii == 1
            value = 'Streams A - B';
        elseif ii == 2 
            value = 'Streams C - D';
        end
    end

    % ************************** SOLVERS **********************************
    function [x, fval, exitflag] = solve_h(k, ii)  % h solver

        % pull previous data and make the guess 
        [h_a_prev, T_w_prev, h_b_prev] = msplit(data, k - 1, ii);
        sol_guess = [h_a_prev, h_b_prev];
        
        % launch solver
        options = optimset('TolX', 1e-7, 'TolFun', 1e-7, ...
		    			   'MaxFunEvals', 1e7, 'MaxIter', 1e7, ...
		    			   'Display', 'iter');
		[x, fval, exitflag] = fsolve(@eqgen, sol_guess, options);

        % compute difference between the equation and zero
	    function eps = eqgen(sol)
            
            % ensure we have column vectors 		
            dhdx_a = zeros(HX_slices, 1);
			dhdx_b = zeros(HX_slices, 1);
            
            % split solution vector 
			h_a = sol(:, 1);
			h_b = sol(:, 2);
            
            % get p and T_w data             
            T_w_a = T_w_prev(:, 1);
            T_w_b = T_w_prev(:, W2);
            p_a = p_data(:, 1, 1, ii);
            p_b = p_data(:, 2, 1, ii);
            
            % dudt
            dudt_a = dudt(h_a, h_a_prev, p_a, fluid{1, ii});
            dudt_b = dudt(h_b, h_b_prev, p_b, fluid{2, ii});
            
            % Q_cond & Q_rad
            Q_cond_aw = Q_cond(h_a, p_a, T_w_a, fluid{1, ii});
            Q_cond_bw = Q_cond(h_b, p_b, T_w_b, fluid{2, ii});
            Q_rad_ae = Q_rad2(h_a, p_a, fluid{1, ii}, ii);
            
			% enthalpy deltas - stream A
			dhdx_a(1) = h_in(1, ii) - h_a(1);
    		dhdx_a(2 : HX_slices) = h_a(1 : HX_slices - 1) - h_a(2 : HX_slices);
            
            % enthalpy deltas - stream B 
            if SWITCH == 1
                dhdx_b(1) = h_in(2, ii) - h_b(1);
                dhdx_b(2 : HX_slices) = h_b(1 : HX_slices - 1) - h_b(2 : HX_slices);
            elseif SWITCH == 0 
                dhdx_b(HX_slices) = h_in(2, ii) - h_b(HX_slices);
                dhdx_b(HX_slices - 1 : -1 : 1) = h_b(HX_slices : -1 : 2)...
                    - h_b(HX_slices - 1 : -1 : 1);
            end 
                       
            % final equations 
			eps_a = - dudt_a * M / t_delta + m * dhdx_a ...
                    + Q_cond_aw + Q_rad_ae;
			eps_b = - dudt_b * M / t_delta + m * dhdx_b ...
                    + Q_cond_bw;
            
			% combine final exit vector
			eps = [eps_a  eps_b];
        end
    end

    function sol = solve_T(i, k, ii)  % T solver, for a fixed HX_slice                          
        
        % get data  
        [~, T_w_prev, ~] = msplit(data, k - 1, ii);
        T_v_prev = T_v_data(:, :, k - 1, ii);
        [h_a, ~, h_b] = msplit(data, k, ii);
        p_a = p_data(:, 1, 1, ii);
        p_b = p_data(:, 2, 1, ii);

        % calculate Q from h 
        Q_cond_aw = Q_cond(h_a(i), p_a(i), T_w_prev(i, 1), fluid{1, ii});   
        Q_cond_bw = Q_cond(h_b(i), p_b(i), T_w_prev(i, W2), fluid{2, ii});       
 
        % launch ode45 
        sol_guess = [T_w_prev(i, :)'; T_v_prev(i, :)'];
        size(sol_guess);
        t_span = [0, t_delta];
        [~, sol] = ode45(@eqgen, t_span, sol_guess);

        % generate equations 
	    function dTdt_w = eqgen(~,T_w)
            
            % ensure we have column vectors 
            dTdx_w = zeros(W3, 1);
            Q_cond = zeros(W3, 1); 
            
            % Calculate K_w and cp
            K_w = K(T_w) * 50;
            cp_w = cp(T_w);
                 
            % Q_cond AT WALL EDGES 
            Q_cond(1) =  K_w(1) / b_x * (T_w(2) - T_w(1)) - Q_cond_aw;
            Q_cond(W2) = K_w(W2) / b_x * (T_w(W2-1) - T_w(W2)) - Q_cond_bw;
            
            % Q_cond AT THE SPLIT
            Q_cond(W2 + 1) = K_w(W2 + 1) / b_x ...
                            * (T_w(W1) + T_w(W2 + 2) - 2 * T_w(W2 + 1));
            Q_cond(W1) = K_w(W1) / b_x ...
                        * (T_w(W2 + 1) + T_w(W1 - 1) + T_w(W1 + 1) - 3 * T_w(W1));
                    
            % Q_cond AT THE END OF THE TAIL 
            Q_cond(W3) = K_w(W3) / b_x * (T_w(W3 - 1) - T_w(W3)); 
                                     
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
            Q_cond(W2 + 2 : W3 - 1) = K_w(W2 + 2 : W3 - 1) / b_x ...
                                      .* dTdx_w(W2 + 2 : W3 - 1); 
            
            % Q_rad & Generate ODEs
            Q_rad = As_w .* F_w * sigma .* (T_ext(1 : W3).^4 - T_w.^4);
            dTdt_w = (Q_rad + Q_cond) ./ (M_w * cp_w);
        end
    end   

    % ************************** PLOTS ************************************
    % STREAM PLOTS
    figure
    subplot(2, 2, 1)
    plot(1 : HX_slices, squeeze(T_data(:, 1, :, 1)), 'r')
    title('Stream A')
    xlabel('Slices')
    ylabel('Temperature')
    
    subplot(2, 2, 2)
    plot(1 : HX_slices, squeeze(T_data(:, N, :, 1)), 'b')
    title('Stream B')
    xlabel('Slices')
    ylabel('Temperature')
    
    subplot(2, 2, 3)
    plot(1 : HX_slices, squeeze(T_data(:, 1, :, 2)), 'r')
    title('Stream C')
    xlabel('Slices')
    ylabel('Temperature')
    
    subplot(2, 2, 4)
    plot(1 : HX_slices, squeeze(T_data(:, N, :, 2)), 'b')
    title('Stream D')
    xlabel('Slices')
    ylabel('Temperature')
    
    % OVERALL PLOTS
    for ii = 1 : 2
        figure
        hold on
        plot(1 : HX_slices, squeeze(T_data(:, 1, :, ii)), 'r')
        plot(1 : HX_slices, squeeze(T_data(:, 2, :, ii)), 'm')
        plot(1 : HX_slices, squeeze(T_data(:, N - 1, :, ii)), 'c')
        plot(1 : HX_slices, squeeze(T_data(:, N, :, ii)), 'b')
        xlabel('Slices')
        ylabel('Temperature')
        name = pull_title(ii);
        title(name);
    end
    
    % WALL PLOTS 
    figure
    subplot(2, 2, 1)
    plot(1 : HX_slices, squeeze(T_data(HX_slices/2, 2 : N - 1 , :, 1)), 'b')
    title('Streams A-B')
    xlabel('Wall Slices')
    ylabel('Temperature')
    
    subplot(2, 2, 2) 
    plot(1 : HX_slices, squeeze(T_data(HX_slices/2, 2 : N - 1, :, 2)), 'b')
    title('Streams A-B')
    xlabel('Wall Slices')
    ylabel('Temperature')
    
    subplot(2, 2, 3)
    plot(1 : W3 - W2, squeeze(T_v_data(HX_slices/2, :, :, 1)), 'g')
    title('Streams A-B')
    xlabel('Wall Slices')
    ylabel('Temperature')
    
    subplot(2, 2, 4)
    plot(1 : W3 - W2, squeeze(T_v_data(HX_slices/2, :, :, 2)), 'g')
    title('Streams C-D')
    xlabel('Wall Slices')
    ylabel('Temperature')
    
end %RN_10_b 
