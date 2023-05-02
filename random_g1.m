function gen_tab = random_g1(Ncand,Nvar,Np_param, P_param, b_param)
% Function that randomly initializes all candidates in the first generation
% gen_tab -- Matrix, [Ncand*Nvar] size  
% Ncand -- number of candidates to genereate
% Nvar -- number of design variables
% Np_param -- [min_teeth, max_teeth]
% P_param -- [min diametral pitch, max diametral pitch]
% b_param -- [min face width, max face width]

% gen_tab = [[Np, P, b]]

% Initialize matrix
gen_tab = zeros(Ncand,Nvar);
szz = size(gen_tab, 1);

% Cycle through each row and fill in with appropriate samples
for k = 1:szz
    % Randomly select starting numbers
    Np_samp = random_sample_from_range(Np_param(1),Np_param(2),1);
    P_samp = random_sample_from_range(P_param(1),P_param(2),1);
    %b_samp = random_sample_from_range(b_param(1),b_param(2),1); %
    %Rememeber, we have to apply constraint on the teeth width also
    
    b_low_lim = 9/P_samp;
    b_high_lim = 14/P_samp;
    %b_samp = random_sample_from_range(b_low_lim,b_high_lim,1);
    b_samp = b_low_lim + ((b_high_lim + b_low_lim) / 16); % Emperical equation, works
    %b_samp = b_low_lim; % Emperical equation, works
    
    % Form a row with these numbers
    rowx = [Np_samp, P_samp, b_samp];
    % Update gen_tab matrix
    gen_tab(k,:) = rowx;
end

end