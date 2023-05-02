function SH = compute_Sh_steel(bhn, CLi, Crr)
    % Function to compute surface endurance strength
    Sfe = (0.4 * bhn) - 10; % ksi, formula for steel from Table 15.5
    Sfe = Sfe * 1000; % conver to psi
    SH = Sfe * CLi * Crr;
end