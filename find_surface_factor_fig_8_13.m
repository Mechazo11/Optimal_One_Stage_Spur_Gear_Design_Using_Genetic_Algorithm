function CS = find_surface_factor_fig_8_13(bhn, surface_finish)
    % Returns the value of surface factor, CS based on Figure 8-13
    
    if (strcmp(surface_finish,"machined"))
        % Using a cubic fit
        % [2.27250393081766e-09,-3.17238039083563e-06,0.000535447551662192,0.776720013477087]
        CS = 2.27250393081766e-09*(bhn^3) -3.17238039083563e-06*(bhn^2) +  0.000535447551662192*bhn + 0.776720013477087;
    end
    
    if (strcmp(surface_finish,"mirror-polished"))
        % Using a cubic fit
        % [2.27250393081766e-09,-3.17238039083563e-06,0.000535447551662192,0.776720013477087]
        CS = 1;
    end
    
    if (strcmp(surface_finish,"commercial-polished"))
        % Using a cubic fit
        CS = -2.532e-09*(bhn^3) + 3.822e-07*(bhn^2)+0.0001364*bhn+0.8815;
    end
end