function [sigma_pinion, sigma_gear] = perform_fatigue_analysis(analysis_type, Np, P, b, gear_ratio, rpm, HP, FOS ,manufac_str, Ko, Km, E_param)
    % Performs the refined gear-tooth-bending strength faituge analysis as
    % per section 15.8
    % gtb - Gear-Tooth-Bending
    % gsd - Geat-Tooth Surface Durability
    
    % Variables
    % Np -- pionion teeth (integer)
    % P -- Diametral Pitch
    % b -- face width
    % rpm -- pinion rpm
    % HP -- horsepower to transmit
    % manufac_str -- string [form_cutting, shaping, high_precision_shaved_ground]
    % needed to determine formula for dynamic velocity factor Kv
    % Ko -- overload factor
    % Km -- mounting factor
    
    Ngg = Np * gear_ratio;
    
    
    if (strcmp(analysis_type, "gtb")==1)
        sigma_pinion = compute_gtb(Np, P, b, rpm, HP, FOS, manufac_str, Ko, Km); % psi
        sigma_gear = compute_gtb(Ngg, P, b, rpm, HP, FOS, manufac_str, Ko, Km); % psi
    end
    
    if (strcmp(analysis_type, "gsd")==1)
        sigma_pinion = compute_gsd_steel_20_deg_pressure(gear_ratio, Np, P, b, rpm, HP, FOS, manufac_str, Ko, Km, E_param); % psi
        sigma_gear = compute_gsd_steel_20_deg_pressure(gear_ratio, Ngg, P, b, rpm, HP, FOS, manufac_str, Ko, Km, E_param); % psi
    end
    
    
end

% -----------------------------------------------------------------------------
function sigma_gtb = compute_gtb(Nteeth, P, b, rpm, HP, FOS ,manufac_str, Ko, Km)
% Local function to compute Gear-Tooth Bending stress for a spur gear
    
    V = (pi * Nteeth * rpm)/(12 * P); % feet per min     
    Ft = (33000 * HP)/ V; % Analytic tangential tooth
    Fd = FOS * Ft; % Design tangential tooth
    
    % Derive dynamic velocity factor Kv
    % HARDCODED FOR THREE CASES
    if(strcmp(manufac_str,"form_cutting"))
        Kv = (600 + V)/600;
    elseif(strcmp(manufac_str,"shaping"))
        %disp("shaping");
        Kv = (1200 + V)/1200;
    else
        %disp("high_precision_shaved_ground");
        Kv = sqrt((78 + sqrt(V))/78);
    end
    
    % Derive J using Np
    if (round(Nteeth)>= 12 & round(Nteeth)<18)
        J = 0.22;
    elseif (round(Nteeth)>= 18 & round(Nteeth)<22)
        J = 0.24;
    else
        J = 0.26;
    end
    
    sigma_gtb = ((Fd * P)/(b * J)) * Kv * Ko * Km; % psi
    
end
% -----------------------------------------------------------------------------

function sigma_H = compute_gsd_steel_20_deg_pressure(gear_ratio, Nteeth, P, b, rpm, HP, FOS ,manufac_str, Ko, Km,E_param)
% Local function to compute Gear-Tooth Surface Durability stress for a spur gear
nu = 0.3; % From Table 15.4a
psi = 20; % degrees
Ep = E_param(1);
Eg = E_param(2);

Cp_deno = ((1 - nu^2)/Ep) + ((1 - nu^2)/Eg);
Cp = sqrt(1 / Cp_deno); %sqrt(psi), a.k.a elastic coefficient
I = ((sin(psi) * cos(psi))/2) * (gear_ratio/(gear_ratio + 1)); % float, called Geometry Factor

V = (pi * Nteeth * rpm)/(12 * P); % feet per min     
Ft = (33000 * HP)/ V; % Analytic tangential tooth load (lb)
Fd = FOS * Ft; % Design tangential tooth (lb)

% Compute Kv
% Derive dynamic velocity factor Kv
    % HARDCODED FOR THREE CASES
if(strcmp(manufac_str,"form_cutting"))
    Kv = (600 + V)/600;
elseif(strcmp(manufac_str,"shaping"))
    %disp("shaping");
    Kv = (1200 + V)/1200;
else
    %disp("high_precision_shaved_ground");
    Kv = sqrt((78 + sqrt(V))/78);
end

d_p = Nteeth / P; % pitch dia of the spur gear

% Compute Hertz stress
sigma_H = Cp * sqrt((Fd/(b * d_p*I))*Kv*Ko*Km);

end
