function [h_data, q_data, T_data, p_data,...
    T_wA_data, T_wB_data, T_Bu_data] = RN_12c
    clc; clear; 
    close all;
    tic
    
    % Variable Data
    HX_UA_in = 30;
    t = 21600; 
    t_delta = 1 * ones(t + 1, 1); 
    em = .2e-04; 
        
    % *********************************************************************
    % ************************** COOL PROP ********************************
    CP_dump = 10; % number of time steps before we dump & re-load the library 
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

    % *************************** HX DATA *********************************
    no = 2; % number of streams, - 
    HX_slices = 20;  % number of slices, -
    HX_length = 10.7;  % cryo-panel length, m 
    slice_per_um = 2; % Number of slices per 10 um of length
    
    % ************************ INITIAL DATA *******************************
    % (dimension j & ii are introduced)
    p_atm = 101325; % define atmospheric pressure to be used as a reference, Pa 
    fluid = 'nitrogen';   
    p_in = 2 * p_atm; % Nitrogen's inlet pressure
    q_in = 0.1; % Nitrogen's inlet quality 
    T_init = 305; % Initital temperature for everything
    h_in = prop_hqp(q_in, p_in, fluid, 'CP'); % convert input to h (nitrogen)
    h_init = h_in;
    
    % mass flow rate for nitrogen 
    v = 100; % L / hr
    v = v / 3600 / 1e3; % convert to m^3 / s
    rho_in = 807; % kg / m^3
    m = rho_in * v;  % kg / s 
    
    % ******************** WALL A1 CALCULATIONS *************************
    L_A1 = HX_length; % m
    H_A1 = 1; % in 
    W_A1 = 0.02; % in 
    % --------------------------------
    H_A1 = H_A1 * 2.54 / 100; % m
    W_A1 = W_A1 * 2.54 / 100; % m
    % --------------------------------
    slice_A1 = round(H_A1 * 1e2 * slice_per_um);
    rho_Cu = 8940; % kg/m^3
    % --------------------------------
    area_A1 = L_A1 * H_A1 / (HX_slices * slice_A1); % m^2
    vol_A1 = area_A1 * W_A1; % m^2
    mass_A1 = vol_A1 * rho_Cu; 
    b_A1 = W_A1 * (L_A1 / HX_slices) / (H_A1 / slice_A1); 
    
    % ******************** WALL A2 CALCULATIONS *************************
    L_A2 = HX_length; % m
    H_A2 = 4.625; % in 
    W_A2 = 0.02; % in 
    % --------------------------------
    H_A2 = H_A2 * 2.54 / 100; % m
    W_A2 = W_A2 * 2.54 / 100; % m
    % --------------------------------
    slice_A2 = round(H_A2 * 1e2 * slice_per_um);
    rho_Cu = 8940; % kg/m^3
    % --------------------------------
    area_A2 = L_A2 * H_A2 / (HX_slices * slice_A2); % m^2
    vol_A2 = area_A2 * W_A2; % m^2
    mass_A2 = vol_A2 * rho_Cu; 
    b_A2 = W_A2 * (L_A2 / HX_slices) / (H_A2 / slice_A2); 

    % ******************** Bulk Head Calculations *************************
    L_Bu = 9; % in 
    H_Bu = ((9 - (1 + 9/16) * 2) * (2 + 1/8) + (1 + 9/16) * 2 * 1) / 9; % in 
    W_Bu = 0.02; % in 
    % --------------------------------
    L_Bu = L_Bu * 2.54 / 100; % m
    H_Bu = H_Bu * 2.54 / 100; % m
    W_Bu = W_Bu * 2.54 / 100; % m
    % --------------------------------
    slice_Bu = round(L_Bu * 1e2 * slice_per_um);
    rho_SS = 7860; % kg/m^3
    % --------------------------------
    area_Bu = L_Bu * H_Bu * 2 * 10 / (HX_slices * slice_Bu); % m^2
    vol_Bu = 1/2 * area_Bu * W_Bu; % m^2
    mass_Bu = vol_Bu * rho_SS; 
    b_Bu = W_Bu * (H_Bu / HX_slices) / (L_Bu / slice_Bu); 
    
    % ******************** Wall Bottom Calculations ***********************
    L_B1 = 10.7; % m 
    H_B1 = 9; % in 
    W_B1 = 0.02; % in 
    % --------------------------------
    H_B1 = H_B1 * 2.54 / 100; % m
    W_B1 = W_B1 * 2.54 / 100; % m
    % --------------------------------
    slice_B1 = round(H_B1 * 1e2 * slice_per_um);
    rho_Cu = 8940; % kg/m^3
    % --------------------------------
    area_B1 = L_B1 * H_B1 / (HX_slices * slice_B1); % m^2
    vol_B1 = area_B1 * W_B1; % m^2
    mass_B1 = vol_B1 * rho_Cu; 
    b_B1 = W_B1 * (L_B1 / HX_slices) / (H_B1 / slice_B1); 

    % Calculate Wall slice intermediates  
    WA1 = slice_A1;
    WA2 = WA1 + slice_A2;
    WB1 = slice_B1;
    WBu = slice_Bu;
    % For Helium
    L_He = [1 + 9/16 - 1/2, 9 - 1 - 9/16 + 1/2] * 2.54 / 100;
    slice_He = round(L_He * 1e2 * slice_per_um); 
    HE1 = slice_He(1);
    HE2 = slice_He(2);
    
    % *************************** STREAM DATA *****************************
    diam = 1 - .02; % from drawing, in
    radius = diam / 2; % in 
    radius = radius * 2.54 / 100; % m
    As = 2 * pi * radius * HX_length; % m^2
    V = pi * radius^2 * HX_length; % m^3
    
    % Scale
    As = As / (HX_slices * (2 * WA2 + WB1 + WBu)); % Divide by total no. wall slices 
    V = V / HX_slices;
    
    % ************************ RADIATION DATA *****************************
    T_ext_init = 300; % External body temperature, K 
    sigma = 5.676e-8; % Stefan-Boltzmann constant, W/m^2 * K^4
    F = .5; % for He-Wall radiation 
        
    % ************************ NOMINAL DATA *******************************
    m_nom = m; 
    T_nom = 100; % nomial temperature, K            
    p_nom = p_atm; % nomial pressure, Pa
    HX_UA_nom_data = HX_UA_in;
    delta_p_nom_data = p_atm/50;
    % calculate other nominal values
    nom_values = nom_calc(HX_UA_nom_data, delta_p_nom_data);
    
    % *************************** PRE-ALLOCATE ****************************
    p_data = zeros(HX_slices, t, no);
    h_data = zeros(HX_slices, t, no);
    T_data = zeros(HX_slices, 2, t, no);
    q_data = zeros(HX_slices, t, no); 
    T_wA_data = zeros(HX_slices, WA2, t, no); 
    T_wB_data = zeros(HX_slices, WB1, t); 
    T_Bu_data = zeros(HX_slices, WBu, t); 
    
    % ************************ BUILT INITIAL DATA *************************
    for ii = 1 : no
        % Grab delta_p_nom
        delta_p_val = nom_values(2);      % 2 = nominal pressure
        
        % BUILD WALL DATA
        As_wA = built_data([area_A1 area_A2], [WA1 WA2 - WA1]);
        M_wA = built_data([mass_A1 mass_A2], [WA1 WA2 - WA1]);
        
        % FIX FIRST TIME STEP 
        % EXTEND INITIAL CONDITIONS TO ALL SLICES 
        T_ext = T_ext_init * ones(max([HX_slices, WA2, WB1, WBu]), 1);
        p_data(:, 1, ii) = p_dist(delta_p_val, p_in);  
        h_data(:, 1, ii) = h_init * ones(HX_slices, 1); 
        T_data(:, 2, 1, ii) = T_init * ones(HX_slices, 1);
        T_wA_data(:, :, 1, ii) = T_init * ones(HX_slices, WA2); 
        T_wB_data(:, :, 1) = T_init * ones(HX_slices, WB1); 
        T_Bu_data(:, :, 1) = T_init * ones(HX_slices, WBu); 
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
            
            % stack
            h_data(:, k, ii) = sol; 
            h_in(2) = sol(end);
            
            % T_wA SOLVER
            for i = 1 : HX_slices
                sol = solve_T_wA(i, k, ii);
                
                %stack in time
                T_wA_data(i, :, k, ii) = sol(end, :); 
            end
            
            % T_He SOLVER
            for i = 1 : HX_slices
                sol = solve_He(i, k, ii);
                
                %stack in time
                T_data(i, 2, k, ii) = sol(end, :); 
            end
        
            % convert h to T (all time steps except k = 1)
            T_data(:, 1, k, ii) = props_thp(h_data(:, k, ii), ...
                p_data(:, 1, ii), fluid);  % T for j = 1
            
            % convert h to q (all time steps except k = 1)
            q_data(:, k, ii) = 1 - props_qhp(h_data(:, k, ii), ...
                p_data(:, 1, ii), fluid);  % T for j = 1
            
            % get T for k = 1
            T_data(:, 1, 1, ii) = props_thp(h_data(:, 1, ii), ...
                p_in * ones(HX_slices, 1), fluid); % T for j = 1
            
            % get q for k = 1
            q_data(:, 1, ii) = 1 - props_qhp(h_data(:, 1, ii), ...
                p_in * ones(HX_slices, 1), fluid); % T for j = 1
        end
        
        % T_wB SOLVER
        for i = 1 : HX_slices
            sol = solve_T_wB(i, k);
                
            % stack in time 
            T_wB_data(i, :, k) = sol(end, :); % wall B 
        end
        
        % T_Bu SOLVER
        for i = 1 : HX_slices
            sol = solve_T_Bu(i, k);
                
            % stack in time 
            T_Bu_data(i, :, k) = sol(end, :); % wall B 
        end
        
        % dump and re-load coolprop
        if rem(k, CP_dump) == 0 
            % check to see if time step is a multiple of CP_dump
            unloadlibrary 'coolprop'
            loadcoolprop; 
        end
    end
    
    % ***************************** H SOLVER ******************************
    function [x, fval, exitflag] = solve_h(k, ii)  % h solver
        
        % pull the current pressure data
        % (dimension k is elimiated) 
        % (dimension ii is elimiated) 
        p = p_data(:, 1, ii);
        
        % pull previous data and make the guess 
        % (dimension k is elimiated) 
        % (dimension ii is elimiated) 
        h_prev = h_data(:, k - 1, ii); 
        
        % launch solver
        options = optimset('TolX', 1e-7, 'TolFun', 1e-7, ...
		    			   'MaxFunEvals', 1e7, 'MaxIter', 1e7, ...
		    			   'Display', 'iter');
		[x, fval, exitflag] = fsolve(@eqgen, h_prev, options);

        % compute difference between the equation and zero
	    function eps = eqgen(h)
            
            T_wA(:, :, 1) = T_wA_data(:, :, k - 1, 1);
            T_wA(:, :, 2) = T_wA_data(:, :, k - 1, 2);
            T_wB = T_wB_data(:, :, k - 1); 
            T_Bu = T_Bu_data(:, :, 1); 
            T_w_prev = [T_wA(:, :, 1), T_wA(:, :, 2), T_wB, T_Bu]; 
            
            % ensure we have column vectors 		
            dhdx = zeros(HX_slices, 1);
            
            % dudt
            dudt = dudt(h, h_prev, p, fluid);
            
            % Q_cond & Q_rad
            Q_cond_wA = Q_cond(h, p, T_wA(:, 1, ii), fluid);   % to top
            if ii == 1
                Q_cond_wB = Q_cond(h, p, T_wB(:, 1), fluid);   % to right 
            elseif ii == 2
                Q_cond_wB = Q_cond(h, p, T_wB(:, end), fluid);   % to left 
            end
            Q_rad_w = Q_rad(h, p, T_w_prev, fluid);
            Q_rad_e = Q_rad(h, p, T_ext(1 : HX_slices), fluid) ...
                * (2 * WA2 + WB1 + WBu);
           
			% enthalpy delta
			dhdx(1) = h_in(ii) - h(1);
    		dhdx(2 : HX_slices) = h(1 : HX_slices - 1) - h(2 : HX_slices);
            
            % calculate mass per cell
            [~,~,rho,~,~] = propsc_LmurhoPr_hp(h, p, fluid);
            M = V * rho;
            
            % final equations 
			eps = - dudt .* M / t_delta(k) + m * dhdx ...
                + Q_cond_wA + Q_cond_wB + Q_rad_w + Q_rad_e;
        end
    end
    
    % ***************************** T_wA SOLVER ***************************
    function sol = solve_T_wA(i, k, ii)  % Tw_A solve, for a fixed HX_slice                           
        
        % get data - T_wA & T_Bu from previous step 
        T_wA_prev = T_wA_data(:, :, k - 1, ii);
        T_Bu_prev = T_Bu_data(:, :, k - 1);
        
        % get data - h & p from the current step 
        p = p_data(:, 1, ii);
        h = h_data(:, k, ii);
        p_l = h_data(:, k, 1); 
        p_r = h_data(:, k, 2);
        h_l = h_data(:, k, 1); 
        h_r = h_data(:, k, 2); 

        % calculate Q from h 
        Q_cond_aw = Q_cond(h(i), p(i), T_wA_prev(i, 1), fluid);   
 
        % launch ode45 
        sol_guess = [T_wA_prev(i, :)'];
        size(sol_guess);
        t_span = [0, t_delta(k)];
        [~, sol] = ode45(@eqgen, t_span, sol_guess);

        % generate equations 
	    function dTdt_wA = eqgen(~, T_wA)
            
            % ensure we have column vectors 
            dTdx_wA = zeros(WA2, 1);
            Q_cond = zeros(WA2, 1); 
            
            % Calculate K_w and cp
            K_wA = K(T_wA, 'Cu');
            cp_wA = cp(T_wA, 'Cu');

            % Q_cond AT STREAM EDGE
            Q_cond(1) =  K_wA(1) * b_A1 * (T_wA(2) - T_wA(1)) - Q_cond_aw;
                    
            % Q_cond AT SHEILD'S EDGE 
            Q_cond(WA2) = K_wA(WA2) * b_A2 * (T_wA(WA2 - 1) - T_wA(WA2));
            
            % Q_cond AT THE SPLIT
            Q_cond(WA1) = K_wA(WA1) * ...
                (b_A1 * (T_wA(WA1 - 1) - T_wA(WA1)) + ...
                b_A2 * (T_wA(WA1 + 1) - T_wA(WA1)) + ...
                b_Bu * (T_Bu_prev(1) - T_wA(WA1))); 
                                     
            % Q_cond IN THE MAIN WALL 
            dTdx_wA(2 : WA1 - 1)        = T_wA(3 : WA1) + T_wA(1 : WA1 - 2) ...
                                        - 2 * T_wA(2 : WA1 - 1);
            dTdx_wA(WA1 + 1 : WA2 - 1)  = T_wA(WA1 + 2 : WA2) + T_wA(WA1 : WA2 - 2) ...
                                        - 2 * T_wA(WA1 + 1 : WA2 - 1);
            Q_cond(2 : WA1 - 1)         = K_wA(2 : WA1 - 1) .* b_A1 ...
                                        .* dTdx_wA(2 : WA1 - 1);
            Q_cond(WA1 + 1 : WA2 - 1)   = K_wA(WA1 + 1 : WA2 - 1) .* b_A2 ...
                                        .* dTdx_wA(WA1 + 1 : WA2 - 1); 
            
            % Q_rad 
            Q_rad_e = As_wA .* F * em * sigma .* (T_ext(1 : WA2).^4 - T_wA.^4);
            Q_rad_l = Q_rad2(h_l(i), p_l(i), T_wA, fluid);
            Q_rad_r = Q_rad2(h_r(i), p_r(i), T_wA, fluid);
            Q_rad = Q_rad_e + Q_rad_l + Q_rad_r;
            
            % Generate ODEs
            dTdt_wA = (Q_rad + Q_cond) ./ (M_wA .* cp_wA);
        end
    end

    % ***************************** T_wB SOLVER ***************************
    function sol = solve_T_wB(i, k)  % Tw_B solver
        % for a fixed HX_slice                          
        
        % get T_wB data from the previous step
        % (dimension k is eliminated)
        % (dimension ii is eliminated)
        T_wB_prev = T_wB_data(:, :, k - 1);
        
        % get h & p data from the current step 
        % (dimension k is elimianted)
        p_l = p_data(:, 1, 1);
        p_r = p_data(:, 1, 2);
        h_l = h_data(:, k - 1, 1);
        h_r = h_data(:, k - 1, 2); 

        % calculate Q from h 
        Q_cond_l = Q_cond(h_l(i), p_l(i), T_wB_prev(i, 1), fluid);   
        Q_cond_r = Q_cond(h_r(i), p_r(i), T_wB_prev(i, WB1), fluid);       
 
        % launch ode45 
        sol_guess = T_wB_prev(i, :)';
        size(sol_guess);
        t_span = [0, t_delta(k)];
        [~, sol] = ode45(@eqgen, t_span, sol_guess);
        
         % generate equations 
	    function dTdt_wB = eqgen(~,T_wB)
            
            % ensure we have column vectors 
            dTdx_wB = zeros(WB1, 1);
            Q_cond = zeros(WB1, 1); 
            
            % Calculate K_w and cp
            K_wB = K(T_wB, 'Cu');
            cp_wB = cp(T_wB, 'Cu');
                 
            % Q_cond AT WALL EDGES 
            Q_cond(1) =  K_wB(1) .* b_B1 * (T_wB(2) - T_wB(1)) - Q_cond_l;
            Q_cond(WB1) = K_wB(WB1) .* b_B1 * (T_wB(WB1 - 1) - T_wB(WB1)) - Q_cond_r; 
                                     
            % Q_cond IN THE MAIN WALL 
            dTdx_wB(2 : WB1 - 1)    = T_wB(3 : WB1) + T_wB(1 : WB1 - 2) ...
                                    - 2 * T_wB(2 : WB1 - 1);
            Q_cond(2 : WB1 - 1)     = K_wB(2 : WB1 - 1) .* b_B1 ...
                                    .* dTdx_wB(2 : WB1 - 1);
            
            % Q_rad 
            Q_rad_e = area_B1 .* F * em * sigma .* (T_ext(1 : WB1).^4 - T_wB.^4);
            Q_rad_l = Q_rad2(h_l(i), p_l(i), T_wB, fluid);
            Q_rad_r = Q_rad2(h_r(i), p_r(i), T_wB, fluid);
            Q_rad = Q_rad_e + Q_rad_l + Q_rad_r; 
            
            % Generate ODEs
            dTdt_wB = (Q_cond + Q_rad) ./ (mass_B1 .* cp_wB);
        end
    end    

    % ***************************** T_Bu SOLVER ***************************
    function sol = solve_T_Bu(i, k)  % Tw_B solver
        % for a fixed HX_slice                          
        
        % get T_wB data from the previous step
        % (dimension k is eliminated)
        % (dimension ii is eliminated)
        T_Bu_prev = T_Bu_data(:, :, k - 1);
        T_wA_prev_l = T_wA_data(:, :, k - 1, 1);
        T_wA_prev_r = T_wA_data(:, :, k - 1, 2);
        T_He_l = T_data(:, 2, k - 1, 1);
        T_He_r = T_data(:, 2, k - 1, 2); 

        % get h & p data from the current step 
        % (dimension k is elimianted)
        p_l = p_data(:, 1, 1);
        p_r = p_data(:, 1, 2);
        h_l = h_data(:, k - 1, 1);
        h_r = h_data(:, k - 1, 2); 
        
        % launch ode45 
        sol_guess = T_Bu_prev(i, :)';
        size(sol_guess);
        t_span = [0, t_delta(k)];
        [~, sol] = ode45(@eqgen, t_span, sol_guess);
        
         % generate equations 
	    function dTdt_Bu = eqgen(~,T_Bu)
            
            % ensure we have column vectors 
            dTdx_Bu = zeros(WBu, 1);
            Q_cond = zeros(WBu, 1); 
            
            % Calculate K_w and cp
            K_Bu = K(T_Bu, 'SS');
            cp_Bu = cp(T_Bu, 'SS');
                 
            % Q_cond AT WALL EDGES 
            Q_cond(1) =  K_Bu(1) .* ...
                 (b_Bu * (T_Bu(2) - T_Bu(1)) + ...
                 b_A1 * (T_wA_prev_l(i, WA1) - T_Bu(1)));
            Q_cond(WBu) = K_Bu(WBu) .* ...
                 (b_Bu * (T_Bu(WBu - 1) - T_Bu(WBu)) + ...
                 b_A1 * (T_wA_prev_r(i, WA1) - T_Bu(WBu))); 
            
            % Q_cond HEAR HE
            Q_cond(HE1) = K_Bu(HE1) * b_Bu ...
                * (T_Bu(HE1 - 1) + T_Bu(HE1 + 1) + T_He_l(i) - 3 * T_Bu(HE1));
            Q_cond(HE2) = K_Bu(HE1) * b_Bu ...
                * (T_Bu(HE2 - 1) + T_Bu(HE2 + 1) + T_He_r(i) - 3 * T_Bu(HE2));
                                     
            % Q_cond IN THE MAIN PART OF BULKHEAD 
            dTdx_Bu(2 : HE1 - 1)        = T_Bu(3 : HE1) + T_Bu(1 : HE1 - 2) ...
                                        - 2 * T_Bu(2 : HE1 - 1);
            dTdx_Bu(HE1 + 1 : HE2 - 1)  = T_Bu(HE1 + 2 : HE2) + T_Bu(HE1 : HE2 - 2) ...
                                        - 2 * T_Bu(HE1 + 1 : HE2 - 1);
            dTdx_Bu(HE2 + 1 : WBu - 1)  = T_Bu(HE2 + 2 : WBu) + T_Bu(HE2 : WBu - 2) ...
                                        - 2 * T_Bu(HE2 + 1 : WBu - 1);
            Q_cond(2 : HE1 - 1)         = K_Bu(2 : HE1 - 1) * b_Bu ...                                    
                                        .* dTdx_Bu(2 : HE1 - 1); 
            Q_cond(HE1 + 1 : HE2 - 1)   = K_Bu(HE1 + 1 : HE2 - 1) * b_Bu ...                                    
                                        .* dTdx_Bu(HE1 + 1 : HE2 - 1);
            Q_cond(HE2 + 1 : WBu - 1)   = K_Bu(HE2 + 1 : WBu - 1) * b_Bu ...
                                        .* dTdx_Bu(HE2 + 1 : WBu - 1);
                            
            % Q_rad 
            Q_rad_e = area_Bu .* F * em * sigma .* (T_ext(1 : WBu).^4 - T_Bu.^4);
            Q_rad_l = Q_rad2(h_l(i), p_l(i), T_Bu, fluid);
            Q_rad_r = Q_rad2(h_r(i), p_r(i), T_Bu, fluid);
            Q_rad = Q_rad_e + Q_rad_l + Q_rad_r;
            
            % Generate ODEs
            dTdt_Bu = (Q_cond + Q_rad) ./ (mass_Bu .* cp_Bu);
        end
    end    

    % ****************************** HE SOLVER ****************************
    function sol = solve_He(i, k, ii)  % HE solver
        % for a fixed HX_slice                          
        
        % get T_wB data from the previous step
        % (dimension k is eliminated)
        % (dimension ii is eliminated)
        T_wA(:, :, 1) = T_wA_data(:, :, k - 1, 1);
        T_wA(:, :, 2) = T_wA_data(:, :, k - 1, 2);
        T_wB = T_wB_data(:, :, k - 1); 
        T_Bu_prev = T_Bu_data(:, :, k - 1);
        T_prev = T_data(:, 2, k - 1, ii);
        T_w_prev = [T_wA(:, :, 1), T_wA(:, :, 2), T_wB, T_Bu_prev];
        
        % launch ode45 
        sol_guess = T_prev(i)';
        t_span = [0, t_delta(k)];
        [~, sol] = ode45(@eqgen, t_span, sol_guess);
        
         % generate equations 
	    function dTdt_He = eqgen(~, T_He)
            
            % Calculate K_w and cp
            h = prop_hqp(T_prev(i), p_atm, 'helium', 'CP');
            [~,~,~,~,cp_He] = propsc_LmurhoPr_hp(h, p_atm, 'helium');
            K_He = 0.142; 
            
            % Mass
            [~,~,rho,~,~] = propsc_LmurhoPr_hp(h, p_atm, 'helium');
            M_He = V * rho;
                 
            % Q_cond AT WALL EDGES 
            Q_cond =  K_He * b_Bu * ...
                (T_Bu_prev(i, slice_He(ii)) - T_He);
            Q_rad = Q_rad3(T_He, T_w_prev(i, :));
            
            % Generate ODEs
            dTdt_He = (Q_cond + Q_rad) ./ (M_He .* cp_He); 
        end
    end    

    % ***************************** FUNCTIONS *****************************
    % FUNCTION THAT CONVERTS DATA INTO PROPER FORMAT (for the wall)
    function value = built_data(data, wall_slices)
        value = zeros(sum(wall_slices), 1); 
        for I = 1 : length(wall_slices)
            value(1 + sum(wall_slices(1 : I - 1)) : sum(wall_slices(1 : I))) = ...
                data(I) * ones(1, wall_slices(I));
        end
    end 
    
    % CALCULATE NOMINAL VALUES
    function values = nom_calc(HX_UA_nom_data, delta_p_nom_data)
      h_nom = prop_htp(T_nom, p_nom, fluid, 'CP');
      [L_nom, mu_nom, rho_nom, Pr_nom] = ...
          propsc_LmurhoPr_hp(h_nom, p_nom, fluid);
      values = [HX_UA_nom_data, delta_p_nom_data, ...
          L_nom, mu_nom, rho_nom, Pr_nom];
    end

    % FUNCTION THAT CALCULATES dudt
    function dudt = dudt(h, h_prev, p, fluid)
		dudt = propsc_uhp(h, p, fluid) - propsc_uhp(h_prev, p, fluid);
    end

    % Calculate HX_UA
    function value = HX_UA(h, p, fluid)
        
        % grab nomimal values 
        HX_UA_nom = nom_values(1);
        L_nom = nom_values(3);
        mu_nom = nom_values(4);
        Pr_nom = nom_values(6);

        % calculate current values 
        [L,mu,~,Pr] = propsc_LmurhoPr_hp(h, p, fluid);

        % scale HX_UA
        value = HX_UA_nom * (L / L_nom) .* (m / m_nom).^.8 .*...
            (mu_nom ./ mu).^.8 .* (Pr / Pr_nom).^(1/3); 
    end

    % FUNCTION THAT CALCULATES Q_cond 
    function Q_cond = Q_cond(h, p, T_w, fluid)
        T = propsc_thp(h, p, fluid);
        T_delta = T_w - T;
        Q_cond = HX_UA(h, p, fluid) * As .* T_delta;
    end

    % FUNCTION THAT CALCULATES Q_rad (stream --> wall) 
    function Q_rad = Q_rad(h, p, T_w, fluid)
        T = propsc_thp(h, p, fluid);
        [~, rad_length] = size(T_w); 
        Q_rad = zeros(HX_slices, rad_length); 
        for I = 1 : rad_length
            Q_rad(:, I) = As * F * em * sigma * (T_w(:, I).^4 - T.^4);
        end
        
        % sum fluxes due to all walls (i.e. across all rows)
        Q_rad = sum(Q_rad, 2);
    end

    % FUNCTION THAT CALCULATES Q_rad (wall --> stream) 
    function Q_rad = Q_rad2(h, p, T_w, fluid)
        [~, rad_length] = size(T_w); 
        T = propsc_thp(h, p, fluid);
        T = T * ones(rad_length, 1); 
        Q_rad = As * F * em * sigma * (T.^4 - T_w.^4);
    end

    % FUNCTION THAT CALCULATES Q_rad (HE --> wall) 
    function Q_rad = Q_rad3(T, T_w)
        [~, rad_length] = size(T_w); 
        Q_rad = zeros(rad_length, 1); 
        for I = 1 : rad_length
            Q_rad(I) = As * F * em * sigma * (T_w(I).^4 - T.^4);
        end
        
        % sum fluxes due to all walls (i.e. across all rows)
        Q_rad = sum(Q_rad); 
    end

    % FUNCTION THAT CALCULATES p_dist 
    function p_dist = p_dist(delta_p, p_in)
        delta_p = delta_p / HX_slices; % delta_p per slice 
        p_dist = p_in : -delta_p : p_in - delta_p * (HX_slices - 1);
    end

    % ******************************* PLOTS *******************************
    figure
    subplot(1, 2, 1)
    plot(1 : HX_slices, q_data(:, :, 1), 'r')
    subplot(1, 2, 2)
    plot(1 : HX_slices, q_data(:, :, 2), 'r')
    
    figure
    subplot(1, 2, 1)
    plot(1 : HX_slices, squeeze(T_data(:, 2, :, 1)), 'r')
    subplot(1, 2, 2)
    plot(1 : HX_slices, squeeze(T_data(:, 2, :, 2)), 'r')
    
    figure
    plot(1 : t + 1, q_data(end, :, end), 'k', 'LineWidth', 2)

    figure
    plot(1 : WA2, squeeze(T_wA_data(HX_slices/2, :, :, 1)), 'b')
    
    figure
    plot(1 : WB1, squeeze(T_wB_data(HX_slices/2, :, :)), 'b')

    figure
    plot(1 : WBu, squeeze(T_Bu_data(HX_slices/2, :, :)), 'b')
end %RN_12c
