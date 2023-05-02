function weight = compute_spur_gear_weight(N, P,b, rho)
    % Function to compute weight of a spur gear
    % Unit - lb
    % N -- Number of teeth (integer)
    % P -- Diametral inch (teeth/in)
    % b -- Face width (in)]
    % rho -- density (lb/in^3)
    weight = (pi/4) * b * rho * ((N/P)^2);
end