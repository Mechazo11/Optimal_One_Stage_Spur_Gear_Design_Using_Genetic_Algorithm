function f_val = obj_f(analysis_type,gen_curr, sn_values, mass_param, rho, ratio, rpm, HP, FOS, manufac_str, Ko, Km, E_param)
% Objective function
% Change this function depending on your optimization problem
% In this project, this is difference between w_max and 

% gen_curr - Ncand * Nvar matrix
% sn_values -- [Sn_pinion, Sn_gear]
% weight_param -- [min_mass, max_mass] % lb
% fval -- [Ncand x 1], column vector of function evaluation

% problem_params = {gear_ratio, pinion_rpm, pow_to_transmit, manufacturing_method, overload_condition, mounting_condition, surface_factor};
% analysis_params = {FOS, kr, kt, kms, CL, CG, CS};

N = size(gen_curr, 1); % Number of rows

% Unpack important variables
m_low = mass_param(1); % lb
m_high = mass_param(2); % lb 
sn_pinion = sn_values(1); % psi
sn_gear = sn_values(2); % psi

% Objective function
f_val = zeros(N,1);

% Compute weight and fitness
for k = 1:N
    % compute_spur_gear_weight(N, P,b, rho)
    Np_in = round(gen_curr(k,1));
    Ng_in = ratio * Np_in; % Derived quantity
    P_in = gen_curr(k,2);
    b_in = gen_curr(k,3);
    
    % Initialize work varialbes
    m_tot = 0.0; 
    m_diff = 0.0;
    
    % Mass of the combined gear pairs
    m_tot = compute_spur_gear_weight(Np_in, P_in, b_in, rho) + ....
    compute_spur_gear_weight(Ng_in, P_in, b_in, rho);
    
    % Guards against impossible mass conditions
    if (m_tot<m_low) 
        m_tot = m_low;
    end
    if (m_tot>m_high) 
        m_tot = m_high;
    end
    
    % Compute first part of fitness function
    m_diff = abs(m_high - m_tot);
    
    % Perform gear-tooth-bending fatigue analysis
    [st_pinion,st_gear] = perform_fatigue_analysis(analysis_type, Np_in, P_in, b_in, ratio, rpm, HP, FOS ,manufac_str, Ko, Km, E_param);
        
    % If design is stable, multiply m_diff by a scalar else multiply it by
    % zero. Zero indicates design does not meet fatigue life constraint
    % Debug message
%     strx = ['sigma_pinion --> ',num2str(sgtb_pinion),' Sn_pinion--> ',num2str(sn_pinion)];
%     disp(strx);
%     strx = ['sigma_gear --> ',num2str(sgtb_gear),' Sn_gear--> ',num2str(sn_gear)];
%     disp(strx);
%     disp(' ');
    
    % Whehter we are using Bending or Surface durability, the logic for
    % both are same here
    if((st_pinion <= sn_pinion) & (st_gear <= sn_gear))
        mm = 1; % Both pinion and gear stress is less than their respective endurance strenghts
    else
        mm = 0.01;
    end
    
    % Final fitness score for this candidate
    m_diff = m_diff * mm;
    
    % Update in f_val vector
    f_val(k,:) = m_diff;
end

end