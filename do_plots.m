function [w_vector, dist_vector] = do_plots(most_fit_matrix, Niter, ratio, rho, material_name)
    % Function to plot weight vs generation after optimization
    w_vector =  zeros(Niter, 1); % Column vector 
    dist_vector =  zeros(Niter, 1); % Column vector 
    
    % Compute weight for each candidate
    for k = 1:Niter
        Npp = most_fit_matrix(k,1); % Number of teeth in pinion
        Ngp = Npp * ratio; % Number of teeth in driven gear
        Pp = most_fit_matrix(k,2); % Diametral pitch
        bp = most_fit_matrix(k,3); % Teeth width
        
        % Weight of the system
        w_vector(k,:) = compute_spur_gear_weight(Npp, Pp, bp, rho) + ....
        compute_spur_gear_weight(Ngp, Pp, bp, rho); % lb  
        %w_vector(k,:) = round(w_vector(k,:),3);
        % Center distance between two gears
        dist_vector(k,:) = ((Npp/Pp) + (Ngp/Pp))/2;
    
    end
    
    % Plots
    figure;
    w_vector = w_vector'; % Make it row vector
    x_vector = linspace(1,Niter, Niter);
    scatter(x_vector, w_vector, 'red', 'filled'); % Plot Weight vs Generation
    ylim([0, max(w_vector)+2]); % Scale Y axis  
    xlabel("Generation (int)");
    ylabel("Weight (lb)");
    tX = ["Weight vs Generation, Material: ",material_name];
    title(tX);
    
%     figure;
%     dist_vector = dist_vector';
%     scatter(x_vector, dist_vector, 'blue', 'filled'); % Plot Weight vs Generation
%     ylim([0, max(dist_vector)+2]); % Scale Y axis
%     xlabel("Generation (int)");
%     ylabel("Center distance (in)");
%     title("Center Distance vs Generation");
    
end