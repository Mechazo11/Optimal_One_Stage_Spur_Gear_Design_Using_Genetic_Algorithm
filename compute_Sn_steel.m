function Sn = compute_Sn_steel(bhn, CL, CG, CS, kr, kt, kms)
    % Script that computes the 10^6 endurance strength based given material
    % properties
    % Works only for steel and needs a Bhn number
    Sn_prime = 0.25 * bhn; % ksi 0.5 * 0.5 * Su
    Sn_prime = Sn_prime * 1000; % conver to psi
    Sn = Sn_prime * CL * CG * CS * kr * kt * kms;
end