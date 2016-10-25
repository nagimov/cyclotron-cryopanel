function [data, T_data, T_v_data, T_wB_data, p_data] = RN_10_d
    % Second version of the 4-stream model 
    % Added flux to "middle walls" 
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
    no = 2; % How many pairs of streams are there
    SWITCH = 0; % 1 for co-current, 0 for counter-current 
    t_delta = .1;  % time step, s
    t = 20;  % number of time steps, -
    HX_slices = 20;  % number of slices, -
    Wall_slices_A = [10, 10, 10];
        % stream A(C) to node, -
        % stream node to stream B(D), -
        % node to tail, -
    Wall_slices_B = 20; % middle sections, -
   
    % CALCULATE WALL SLICE INTERMEDIATES 
    WA1 = Wall_slices_A(1); % from stream A(C) to node
    WA2 = WA1 + Wall_slices_A(2); % from node to stream B(D)
    WA3 = WA2 + Wall_slices_A(3); % from node to tail  
    WB = Wall_slices_B; % in-between streams 
    N = 1 + WA2 + 1; 
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
    T_in = [150, 100;               % stream A, stream C; temperature, K
            200, 150];             % stream B, stream D; temperature, K
    h_in = hcalc(T_in, p_in, fluid); % convert input to h     
    T_wA_init = [175, 125]; % initial wall temperature (left, right), K, -   
    T_wB_init = [125, 175]; % initial wall temperature (top, bottom), K, - 
    T_v_init = [175, 125]; % initial wall tail temperature (left, right), K, - 
    
    % ************************ RADIATION DATA *****************************
    T_ext_init = 300; % External body temperature, K 
    sigma = 5.676e-8; % Stefan-Boltzmann constant, W/m^2 * K^4 
    As_wA_data = [1, 1, 10]; % Surface area, m^2 
    F_wA_data = [1, 1, 1];  % Radiation view factor, -
        % stream A(C) to node, 
        % stream node to stream B(D), 
        % node to tail    
    As_wB_data = [1, 1]; % Surface area, m^2 
    F_wB_data = [1, 1];  % Radiation view factor, -    
        % top, bottom 

    % ************************ NOMINAL DATA *******************************
    m_nom = 1; % nominal mass flow rate, kg/s
    T_nom = 100; % nomial temperature, K            
    p_nom = p_atm; % nomial pressure, Pa
    HX_UA_nom_data = {'nitrogen', 2500; ...
                        'helium', 2500}; % nomial HX coefficient, W/K
    delta_p_nom_data = {'helium', p_atm / 3; ... 
                        'nitrogen', p_atm / 3}; % nomial pressure drop, Pa
    nom_values = {'helium', nom_calc('helium'); ...
                   'nitrogen', nom_calc('nitrogen')}; % calculate other nominal values
               
    % ************************ PRE-ALLOCATE *******************************
    % The data matrices keep the whole solution
    % i are the HX_slices
    % j are the Wall_slices_A
    % k are the time steps 
    % ii is the side of the cryo-panel 
    data = zeros(HX_slices, N, t + 1, no);
    
    % data matrix converted into T
    T_data = zeros(HX_slices, N, t + 1, no);
    
    % stores wall tail T
    T_v_data = zeros(HX_slices, WA3-WA2, t + 1, no);    
    % stores wall B data 
    T_wB_data = zeros(HX_slices, WB, t + 1, no); 
    
    % stores pressure 
    p_data = zeros(HX_slices, 2, t + 1, no);
    h_0 = zeros(HX_slices, 2, no);
    p_0 = zeros(HX_slices, 2, no);

    % ************************ BUILT INITIAL DATA *************************
    for ii = 1 : no
    
        % Built rad data for walls   
        As_wA = built_rad_data(As_wA_data, Wall_slices_A);
        F_wA = built_rad_data(F_wA_data, Wall_slices_A);
        As_wB = built_rad_data(As_wB_data, Wall_slices_B);
        F_wB = built_rad_data(F_wB_data, Wall_slices_B); 
        
        % EXTEND INITIAL CONDITIONS TO ALL SLICES 
        % (dimension i is introduced)
        % (dimension j is extended)
        h_0(:, 1, ii) = h_in(1, ii) * ones(HX_slices, 1); 
        h_0(:, 2, ii) = h_in(2, ii) * ones(HX_slices, 1);
        p_0(:, 1, ii) = p_in(1, ii) * ones(HX_slices, 1);
        p_0(:, 2, ii) = p_in(2, ii) * ones(HX_slices, 1);
        T_wA_0 = T_wA_init(ii) * ones(HX_slices, WA2);
        T_wB_0(:, :, 1) = T_wB_init(1) * ones(HX_slices, WB);
        T_wB_0(:, :, 2) = T_wB_init(2) * ones(HX_slices, WB);
        T_v_0 = T_v_init(ii) * ones(HX_slices, WA3 - WA2);
        T_ext = T_ext_init * ones(max(HX_slices, WA3), 1); 
                       
        % COMBINE INITIAL CONDITIONS 
        % (dimension k is introduced)
        data(:, :, 1, ii) = [h_0(:, 1, ii), T_wA_0, h_0(:, 2, ii)];
        p_data(:, :, 1, ii) = p_0(:, : , ii); 
        T_v_data(:, :, 1, ii) = T_v_0; 
        T_wB_data(:, :, 1, 1) = T_wB_0(:, :, 1);
        T_wB_data(:, :, 1, 2) = T_wB_0(:, :, 2);
        
    end 
        
    % ************************** LAUNCH SOLVERS ***************************
    for k = 2 : t + 1
        for ii = 1 : no
            
            % h SOLVER
            [sol, fval, exitflag] = solve_h(k, ii);
            disp(['Time step completed ' num2str(k - 1) '(' num2str(ii) ')'...
                ' Time ' num2str(toc/60) ' min '...
                ' Fval ' num2str(sum(fval)) ...
                ' Exit Flag ' num2str(exitflag)])
            
            % stack in time 
            data(:, 1, k, ii) = sol(:, 1); 
            data(:, N, k, ii) = sol(:, 2); 
            
            % T_wA & T_wB SOLVER
            for i = 1 : HX_slices
                solA = solve_T_wA(i, k, ii);
                solB = solve_T_wB(i, k, ii);
                
                %stack in time
                data(i, 2 : N - 1, k, ii) = solA(end, 1 : WA2); % main wall
                T_v_data(i, :, k, ii) = solA(end, WA2 + 1 : WA3); % wall tail
                T_wB_data(i, :, k, ii) = solB(end, :); % wall B 
            end
            
            % convert h to T (all time steps except k = 1)
            T_data(:, 1, k, ii) = props_thp(data(:, 1, k, ii), ...
                p_data(:, 1, 1, ii), fluid{1, ii});  % T for j = 1
            T_data(:, N, k, ii) = props_thp(data(:, N, k, ii), ...
                p_data(:, 1, 1, ii), fluid{2, ii});  % T for j = N
            
            % get T for k = 1
            T_data(:, 1, 1, ii) = props_thp(data(:, 1, 1, ii), ...
                p_data(:, 1, 1, ii), fluid{1, ii}); % T for j = 1
            T_data(:, N, 1, ii) = props_thp(data(:, N, 1, ii), ...
                p_data(:, 1, 1, ii), fluid{2, ii});  % T for j = N
            
            % dump and re-load coolprop
            if rem(k, CP_dump) == 0 
                % check to see if time step is a multiple of CP_dump
                unloadlibrary 'coolprop'
                loadcoolprop; 
            end 
        end
    end
    
    % copy wall temperatures to T_data matrix
    T_data(:, 2 : N - 1, :, :) = data(:, 2 : N - 1, :, :); % for j NOT
    plots
    
    % ***************************** FUNCTIONS *****************************
    % Enthalpy conversion  
    function h = hcalc(T, p, fluid)
        h = zeros(2); 
        for j = 1 : 2
            for ii = 1 : no
                if T(j, ii) <= 1 % this means we have q as an input
                    h(j, ii) = prop_hqp(T(j, ii), p(j, ii), fluid{j, ii}, 'CP');
                else % this means we have T as an input 
                    h(j, ii) = props_htp(T(j, ii), p(j, ii), fluid{j, ii});
                end
            end
        end 
    end 

    % Split data matrix 
    function [h_a, T_wA, h_b] = msplit(matrix, k, ii)
        h_a = matrix(:, 1, k, ii);
        T_wA = matrix(:, 2 : N - 1, k, ii);
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
    function value = built_rad_data(rad_data_w, wall_slices)
        value = zeros(sum(wall_slices), 1); 
        for I = 1 : length(wall_slices)
            value(1 + sum(wall_slices(1 : I - 1)) : sum(wall_slices(1 : I))) = ...
                rad_data_w(I) * ones(1, wall_slices(I));
        end
    end

    % FUNCTION THAT CALCULATES Q_cond 
    function Q_cond = Q_cond(h, p, T_wA, fluid)
        T = propsc_thp(h, p, fluid);
        T_delta = T_wA - T;
        Q_cond = HX_UA(h, p, fluid) / HX_slices .* T_delta;
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
        [h_a_prev, T_wA_prev, h_b_prev] = msplit(data, k - 1, ii);
        T_wB_prev = T_wB_data(:, :, k - 1, :); 
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
            
            % get p, T_wA and T_wB data             
            T_wA_a = T_wA_prev(:, 1);
            T_wA_b = T_wA_prev(:, WA2);
            if ii == 1
                T_wB_a = T_wB_prev(:, 1, 1);
                T_wB_b = T_wB_prev(:, 1, 2);
            elseif ii == 2
                T_wB_a = T_wB_prev(:, WB, 1);
                T_wB_b = T_wB_prev(:, WB, 2);
            end
            p_a = p_data(:, 1, 1, ii);
            p_b = p_data(:, 2, 1, ii);
            
            % dudt
            dudt_a = dudt(h_a, h_a_prev, p_a, fluid{1, ii});
            dudt_b = dudt(h_b, h_b_prev, p_b, fluid{2, ii});
            
            % Q_cond & Q_rad
            Q_cond_awA = Q_cond(h_a, p_a, T_wA_a, fluid{1, ii});
            Q_cond_bwA = Q_cond(h_b, p_b, T_wA_b, fluid{2, ii});
            Q_cond_awB = Q_cond(h_a, p_a, T_wB_a, fluid{1, ii});
            Q_cond_bwB = Q_cond(h_b, p_b, T_wB_b, fluid{2, ii});
            
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
                    + Q_cond_awA + Q_cond_awB; 
			eps_b = - dudt_b * M / t_delta + m * dhdx_b ...
                    + Q_cond_bwA + Q_cond_bwB;
            
			% combine final exit vector
			eps = [eps_a  eps_b];
        end
    end

    function sol = solve_T_wA(i, k, ii)  % Tw_A solver, for a fixed HX_slice                          
        
        % get data  
        [~, T_wA_prev, ~] = msplit(data, k - 1, ii);
        T_v_prev = T_v_data(:, :, k - 1, ii);
        [h_a, ~, h_b] = msplit(data, k, ii);
        p_a = p_data(:, 1, 1, ii);
        p_b = p_data(:, 2, 1, ii);

        % calculate Q from h 
        Q_cond_aw = Q_cond(h_a(i), p_a(i), T_wA_prev(i, 1), fluid{1, ii});   
        Q_cond_bw = Q_cond(h_b(i), p_b(i), T_wA_prev(i, WA2), fluid{2, ii});       
 
        % launch ode45 
        sol_guess = [T_wA_prev(i, :)'; T_v_prev(i, :)'];
        size(sol_guess);
        t_span = [0, t_delta];
        [~, sol] = ode45(@eqgen, t_span, sol_guess);

        % generate equations 
	    function dTdt_wA = eqgen(~,T_wA)
            
            % ensure we have column vectors 
            dTdx_wA = zeros(WA3, 1);
            Q_cond = zeros(WA3, 1); 
            
            % Calculate K_w and cp
            K_wA = K(T_wA) * 50;
            cp_wA = cp(T_wA);
                 
            % Q_cond AT WALL EDGES 
            Q_cond(1) =  K_wA(1) / b_x * (T_wA(2) - T_wA(1)) - Q_cond_aw;
            Q_cond(WA2) = K_wA(WA2) / b_x * (T_wA(WA2 - 1) - T_wA(WA2)) - Q_cond_bw;
            
            % Q_cond AT THE SPLIT
            Q_cond(WA2 + 1) = K_wA(WA2 + 1) / b_x ...
                            * (T_wA(WA1) + T_wA(WA2 + 2) - 2 * T_wA(WA2 + 1));
            Q_cond(WA1) = K_wA(WA1) / b_x ...
                        * (T_wA(WA2 + 1) + T_wA(WA1 - 1) + T_wA(WA1 + 1) - 3 * T_wA(WA1));
                    
            % Q_cond AT THE END OF THE TAIL 
            Q_cond(WA3) = K_wA(WA3) / b_x * (T_wA(WA3 - 1) - T_wA(WA3)); 
                                     
            % Q_cond IN THE MAIN WALL 
            dTdx_wA(2 : WA1 - 1)      = T_wA(3 : WA1) + T_wA(1 : WA1 - 2) ...
                                      - 2 * T_wA(2 : WA1 - 1);
            dTdx_wA(WA1 + 1 : WA2 - 1) = T_wA(WA1 + 2 : WA2) + T_wA(WA1 : WA2 - 2) ...
                                       - 2 * T_wA(WA1 + 1 : WA2 - 1);
            Q_cond(2 : WA1 - 1)      = K_wA(2 : WA1 - 1) / b_x ...
                                     .* dTdx_wA(2 : WA1 - 1);
            Q_cond(WA1 + 1 : WA2 - 1) = K_wA(WA1 + 1 : WA2 - 1) / b_x ...
                                      .* dTdx_wA(WA1 + 1 : WA2 - 1); 
            
            % Q_cond IN THE TAIL
            dTdx_wA(WA2 + 2 : WA3 - 1) = T_wA(WA2 + 3 : WA3) + T_wA(WA2 + 1 : WA3 - 2) ...
                                       - 2 * T_wA(WA2 + 2 : WA3 - 1);
            Q_cond(WA2 + 2 : WA3 - 1) = K_wA(WA2 + 2 : WA3 - 1) / b_x ...
                                      .* dTdx_wA(WA2 + 2 : WA3 - 1); 
            
            % Q_rad & Generate ODEs
            Q_rad = As_wA .* F_wA * sigma .* (T_ext(1 : WA3).^4 - T_wA.^4);
            dTdt_wA = (Q_rad + Q_cond) ./ (M_w * cp_wA);
        end
    end

    function sol = solve_T_wB(i, k, ii)  % Tw_B solver, for a fixed HX_slice                          
        
        % get data  
        T_wB_prev = T_wB_data(:, :, k - 1, ii);
        if ii == 1
            [h_l, ~, ~] = msplit(data, k, 1);
            [h_r, ~, ~] = msplit(data, k, 2);
            p_l = p_data(:, 1, 1, 1);
            p_r = p_data(:, 1, 1, 2);
        elseif ii == 2
            [~, ~, h_l] = msplit(data, k, 1);
            [~, ~, h_r] = msplit(data, k, 2);
            p_l = p_data(:, 2, 1, 1);
            p_r = p_data(:, 2, 1, 2);
        end 
        
        % calculate Q from h 
        Q_cond_l = Q_cond(h_l(i), p_l(i), T_wB_prev(i, 1), fluid{ii, 1});   
        Q_cond_r = Q_cond(h_r(i), p_r(i), T_wB_prev(i, WB), fluid{ii, 2});       
 
        % launch ode45 
        sol_guess = T_wB_prev(i, :)';
        size(sol_guess);
        t_span = [0, t_delta];
        [~, sol] = ode45(@eqgen, t_span, sol_guess);
        
         % generate equations 
	    function dTdt_wB = eqgen(~,T_wB)
            
            % ensure we have column vectors 
            dTdx_wB = zeros(WB, 1);
            Q_cond = zeros(WB, 1); 
            
            % Calculate K_w and cp
            K_wB = K(T_wB) * 50;
            cp_wB = cp(T_wB);
                 
            % Q_cond AT WALL EDGES 
            Q_cond(1) =  K_wB(1) / b_x * (T_wB(2) - T_wB(1)) - Q_cond_l;
            Q_cond(WB) = K_wB(WB) / b_x * (T_wB(WB - 1) - T_wB(WB)) - Q_cond_r; 
                                     
            % Q_cond IN THE MAIN WALL 
            dTdx_wB(2 : WB - 1) = T_wB(3 : WB) + T_wB(1 : WB - 2) ...
                                - 2 * T_wB(2 : WB - 1);
            Q_cond(2 : WB - 1)  = K_wB(2 : WB - 1) / b_x ...
                                .* dTdx_wB(2 : WB - 1);
            
            % Q_rad & Generate ODEs
            Q_rad = As_wB .* F_wB * sigma .* (T_ext(1 : WB).^4 - T_wB.^4);
            dTdt_wB = (Q_cond + Q_rad) ./ (M_w * cp_wB);
        end
    end    
        
    % ************************** PLOTS ************************************
    function plots
        
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
    
    % WALL PLOTS 
    FigHandle = figure;
    set(FigHandle, 'Position', [500, 50, 600, 800]);
    subplot(3, 2, 1)
    plot(1 : HX_slices, squeeze(T_data(HX_slices/2, 2 : N - 1 , :, 1)), 'b')
    title('Streams A-B')
    xlabel('Wall Slices')
    ylabel('Temperature')
    
    subplot(3, 2, 2) 
    plot(1 : HX_slices, squeeze(T_data(HX_slices/2, 2 : N - 1, :, 2)), 'b')
    title('Streams C-D')
    xlabel('Wall Slices')
    ylabel('Temperature')    
    
    subplot(3, 2, 3)
    plot(1 : WB, squeeze(T_wB_data(HX_slices/2, :, :, 1)), 'b')
    title('Streams A-C')
    xlabel('Wall Slices')
    ylabel('Temperature')
    
    subplot(3, 2, 4)
    plot(1 : WB, squeeze(T_wB_data(HX_slices/2, :, :, 2)), 'b')
    title('Streams B-D')
    xlabel('Wall Slices')
    ylabel('Temperature')
    
    subplot(3, 2, 5)
    plot(1 : WA3 - WA2, squeeze(T_v_data(HX_slices/2, :, :, 1)), 'g')
    title('Streams A-B, Tail')
    xlabel('Wall Slices')
    ylabel('Temperature')
    
    subplot(3, 2, 6)
    plot(1 : WA3 - WA2, squeeze(T_v_data(HX_slices/2, :, :, 2)), 'g')
    title('Streams C-D, Tail')
    xlabel('Wall Slices')
    ylabel('Temperature')
    
    % 3D PLOT PREP
    [X, Y] = meshgrid(1 : HX_slices, 1 : N); % get a mesh grid
    time = [2, round(t / 2), round(3 * t / 4), t + 1]; % which times to plot
    colormap autumn 
    
    for ii = 1 : no

        % OVERALL PLOTS
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
        
        % 3D PLOT
        FigHandle = figure;
        set(FigHandle, 'Position', [150, 50, 1200, 700]);
        for I = 1 : 4 %4 sub-plots
            T_max = max(T_in(:, ii));
            T_min = min(T_in(:, ii));
            subplot(2, 2, I)
            surf(X, Y, T_data(:, :, time(I), ii)')
            xlabel('Length Slices')
            ylabel('Wall Slices')
            zlabel('Temperature')
            title(['Time step ' num2str(time(I))])
            axis([1 HX_slices 1 N T_min T_max])
            shading interp;
            caxis manual
            caxis([T_min T_max]); % fix color distribution to always be between Tmin, Tmax
            colorbar 
        end 
    end
    end
end %RN_10_d
