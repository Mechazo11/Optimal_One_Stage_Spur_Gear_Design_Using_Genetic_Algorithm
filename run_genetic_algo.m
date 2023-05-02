% --------------------------------------------------------------------------
function [gen_future,most_fit_g1] = run_genetic_algo(analysis_type, gen_count, gen_current, fatigue_lifes ,material_params, problem_params, analysis_params, opti_params ,bit_count, Nvar, Ncand, rho)
    % Master function, heavily modified for 3 variable problem
    % Return new generaiton and most fit candidate of this iteration
    
    % Necessary values used throughout the code
    % analysis_type = {gtb, gsd} % gtb - Gear tooth bending strength, gsd -
    % Gear tooth surface durability
    m_min = opti_params{1,7};
    m_max = opti_params{1,8};
    Sn_pinion = fatigue_lifes{1,1};
    Sn_gear = fatigue_lifes{1,2};
    rho = material_params{1,1};
    E_pinion = material_params{1,2};
    E_gear = material_params{1,3};
    ratio = problem_params{1,1};
    rpm = problem_params{1,2};
    HP =  problem_params{1,3};
    manufac_str = problem_params{1,4};
    Ko = problem_params{1,5};
    Km = problem_params{1,6};
    FOS = analysis_params{1,1};
    var_lower_limits = [opti_params{1,1}, opti_params{1,3}, opti_params{1,5}]; % min_Np, min_DP, min_b
    var_upper_limits = [opti_params{1,2}, opti_params{1,4}, opti_params{1,6}]; % max_Np, max_DP, max_b
    
    % Step 0: Check if we have more than Ncand rows in current generation.
    % Bring input list back to Ncand number of rows
       
    % Count number of rows (activate this later)
    row_in = size(gen_current, 1);
    if (row_in > Ncand)
        [gen_current] = bring_back_N(gen_current, Nvar, [Sn_pinion, Sn_gear], [m_min, m_max], rho, ratio, rpm, HP, FOS, manufac_str, Ko, Km);
    end
    
    % Step 0.1
    num_to_gen = Ncand * Nvar;
    rand_ls = random_generator(num_to_gen, 0, 0.98); % Generate prob between [0, 0.98]
    
    % Step 1: Evalute objective function and find cumulative fit
    % Here either Bending or Surface Durability test is performed
    [fit_g1, cumu_fit_g1] = eval_obj(analysis_type, gen_current, [Sn_pinion, Sn_gear],.... 
    [m_min, m_max], rho, ratio, rpm, HP, FOS,.... 
    manufac_str, Ko, Km, [E_pinion, E_gear]);
    
    % Step 2: Find the most fit individual and remember it (Elitisim -- Simply copy and paste)
    % By remembering, we record the values of Np, P and b for this
    % candidate
    fit_mat = [gen_current, transpose(fit_g1)];
    most_fit_g1 = find_most_fit(fit_mat);
    
    % Step 3 Evaluate fraction of fitness and cumulative probability
    [frac_fit_g1, cumu_prob_g1] = eval_fraction_fitness(fit_g1, cumu_fit_g1);
    
 
%     %----------------------------------------------------------------
%     % Step 4 Select N random numbers from pool for genetic cross over
%     % datasample -- MATLAB's builtin function
      comp_rand = datasample(rand_ls, Ncand,'Replace',false);
%     
%     % Step 5 Determine mate using wheel of fortune
%     %mate_matrix = find_mates_without_replacement(X1, cumu_prob_g1,
%     comp_rand); % No repetation of mates is allowed
      mate_matrix = find_mates_with_replacement(gen_current, cumu_prob_g1, comp_rand); % Allow a candiate to be repeated in the mating pool
%     
%     % Step 6 Reshape to have candiates form couple per row
      mmx = reshape_long_row(mate_matrix);
%     
%     % Step 7 Perform cross over
%     % How many couples do we have
      mmx_couple_count = size(mmx,1);
      
%     % Create (N/2) random numbers, row vector
      cross_rand = random_generator(mmx_couple_count, 0, 0.98);
%     
%     % Perform cross over
      % HARDCODED, problem specific
      % Constraint on gear tooth width is applied here
      gen_future = create_new_gen_spur_gears(mmx,bit_count,var_upper_limits,var_lower_limits,Nvar);
%     
%     % Step 8 Reshape new generation matrix to N by n_feature candiates
      gen_future = reshape_to_Ncand_by_Nvar(gen_future, Ncand, Nvar);
%     
%     % Step 9 Perform elitism, copy-paste the strongest guy from the last generation
      %gen_future = [gen_future;most_fit_g1];
%     
end

% ---------------------------------------------------------------------------
function [X1_best] = bring_back_N(gen_current, Nvar, sn_values, mass_param, rho, ratio, rpm, HP, FOS, manufac_str, Ko, Km)
    % HARDCODED FOR THIS PROBLEM
    % Upgraded for 3 variable problem
    % Calculate fitness and take out the weakest member
    [fit_g1, cumu_fit_g1] = eval_obj(gen_current, sn_values, mass_param, rho, ratio, rpm, HP, FOS, manufac_str, Ko, Km);
    fit_g1 = transpose(fit_g1); % Turn row into column vector
    X3 = [gen_current, fit_g1]; % (N by 3)
    [M,I] = min(X3(:,3));
    X3_1 = X3(1:(I-1),:);
    X3_2 = X3((I+1):end,:);
    
    % Drop (n_feature + 1) column we will add it back later
    X3_1 = X3_1(:,(1:Nvar));
    X3_2 = X3_2(:,(1:Nvar));
    X3 = [X3_1;X3_2];
    
    % Return candidates with N rows
    X1_best = X3;
    
end
% -----------------------------------------------------------------------------

% ---------------------------------------------------------------------------
function [fitness, cumu_fit]  = eval_obj(analysis_type, gen_current, sn_values, mass_param, rho, ratio, rpm, HP, FOS, manufac_str, Ko, Km, E_params) 
    % Evaluates objective function and adds up all fitness values
    
    % Evaluate fitness of each candidate in this generation, output as a
    % column vector
    fitness = obj_f(analysis_type, gen_current, sn_values, mass_param, rho, ratio, rpm, HP, FOS, manufac_str, Ko, Km, E_params);
    fitness = transpose(fitness); % Converting to a row vector
    
    % Add up all fitness values
    cumu_fit = sum(fitness,2); % Value 2 tells MATLAB to sum along all the rows
end
% -----------------------------------------------------------------------------

% ---------------------------------------------------------------------------
function [most_fit] = find_most_fit(fit_mat) 
    B = sortrows(fit_mat,3, 'descend'); % Most fit individual on top
    most_fit = B(1,[1:3]); % The three parameters of the best fit
end
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
function [frac_to_fit, cumu_prob] = eval_fraction_fitness(fitness, cumu_fit)
    frac_to_fit = fitness./cumu_fit;
    cumu_prob = cumsum(frac_to_fit,2);
end
% ---------------------------------------------------------------------------

% ---------------------------------------------------------------------------
function X22 = find_mates_with_replacement(gen_current, cumu_prob, comp_rand)
    % Allow a candiate to be repeated in the mating pool
    % Note, both cumu_prob and comp_rand are row vectors
    % They are needed to be flipped to column vector
    % Changelog, fixed issue regarding choosing the column containing
    % cumulative probability and 
    
    % Debug
    %comp_rand = [0.9117, 0.5981, 0.5866, 0.8035, 0.5276, 0.3441];
    
    X11 = [gen_current, (cumu_prob).', (comp_rand).']; % Big matrix
    NX = size(X11,1); % Number of rows i.e number of candiates
    Nvar = size(gen_current,2);
    X22 = zeros(NX,Nvar); % Initialize?
    num_to_iter = size(X11,1);

    % Scratch pad variable
    row_prob = 0;
    rand_prob = 0;
    test_idx = 1;
    cnt = 1;
    
    %X11 % Debug
    
    % Main loop
    for ii = 1:num_to_iter
        rand_prob = X11(ii,end); % The last row contains the random probabilties
        %fprintf("Rand prob now --> %f\n",rand_prob);
        
        % Look through each cumulative probability
        for kk = 1: num_to_iter
            % Load each cumulative probability sequentially
            row_prob = X11(kk,end-1);
            %fprintf("Cumulative prob now --> %f\n",row_prob);
            
            % Check if current cumulative prob is bigger than current row's
            % random probability
            if (row_prob > rand_prob)
                X22(ii,:) = gen_current(kk,:);
                % Reset for next random number
                rand_prob = 0;
                row_prob = 0;
                cnt = 1;
                break; % We got a hit
            else
                cnt = cnt + 1; % Dummy to prevent error
                if (cnt<num_to_iter)
                    continue % Go to next random probability value
                else
                    % Guard against overflow
                    % change random probability value
                    % Exceeded num_iter count
                    % Reset cnt
                    % reset kk
                    % continue again
                    rand_prob = 0.111;
                    cnt = 1;
                    kk = 1;
                    continue
                end
            end
        end
    end
    %X22 % What is X22 now?
end
% ---------------------------------------------------------------------------

% ---------------------------------------------------------------------------
function mxx_reshape = reshape_long_row(mate_matrix)
% Function which reshapes Ncand by Nvar to (Ncand/2) by (Nvar*2) matrix
mmx = mate_matrix;
mmx = mmx'; % [Nvar x Ncand]
mmx = reshape(mmx,1,[]); % All elements flattened into a row vector
num_count = size(mmx,2);
szz1 = size(mate_matrix,1);
szz2 = size(mate_matrix,2);
ppx = zeros((szz1/2),(szz2 * 2)); % Temporary array to hold (N/2) by (n_feature*2)

% Record number of rows and columns of the reshape matrix
ppx_row = size(ppx,1);
ppx_col = size(ppx,2);
ctn = 1;
  % This nested loop is actually extracting the genes's for each candidate
  % from the flattened matrix
  for pp_row = 1: ppx_row % Row selector
      for pp_col = 1:ppx_col % column selector
          ppx(pp_row, pp_col) = mmx(1,ctn); 
          ctn = ctn + 1;
       end
  end
mxx_reshape = ppx;
end
% ---------------------------------------------------------------------------

% -----------------------------------------------------------------------------
function X22 = find_mates_without_replacement(X_current, cumu_prob, comp_rand)
    % Does not allow a candidate to be repeated in the mating pool
    % TODO does not work with this example 03/29/23
    
    X11 = [X_current, (cumu_prob).', (comp_rand).']; % Big matrix
    NX = size(X11,1); % Number of rows i.e number of candiates
    X22 = zeros(NX,2);
    num_to_iter = size(X11,1);

    % Scratch pad variable
    row_prob = 0;
    rand_prob = 0;
    test_idx = 1;
    cnt = 1;

    % Counter to keep track of how many counts we did
    for ii = 1: num_to_iter
        row_prob = X11(ii,3); % Choose the cumu_probability for this candidate
        %fprintf("\nCumu prob now --> %f\n",row_prob);
        %fprintf("------------------------\n");
        
        % Look thru each "random" probability sequentially
        for kk = 1: num_to_iter
            rand_prob = X11(kk,4); % Get probability for position at kk
            %fprintf("Rand prob now --> %f\n",rand_prob);
            %%%Test if current cumu probability is greater than current value of random probability
            if (row_prob >= rand_prob)
                X22(kk,:) = X_current(ii,:);
                X11(kk,4) = 1.1; % This will prevent this position random number to evaluated again in next round
                
                %fprintf("-----------------------------------------------------------\n");
                %fprintf("Assigned G_%d position %d going to next candidate\n", ii,kk);
                %fprintf("-----------------------------------------------------------\n");
                rand_prob = 0;
                row_prob = 0;
                cnt = 1;
                break % We stop the inner loop and move onto next candidate
            else
                cnt = cnt + 1; % Dummy to prevent error
                if (cnt<num_to_iter)
                    continue % Go to next random probability value
                else
                    % change actual probability value of this candidate to 1
                    % Exceeded num_iter count
                    % Reset cnt
                    % reset kk
                    % continue again
                    row_prob = 0.99999;
                    cnt = 1;
                    kk = 1;
                    continue
                end
            end
        end
    end
% X22 now is in N x n_feature form. We will reshape them into (N/2,
% n_feature *2) for easier computation of next generation design candidates
end
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
function [mmx_reshape] = reshape_to_Ncand_by_Nvar(new_gen, N, n_feature)
    % Reshape new_gen matrix to 6 by 2
    tempo_mat = zeros(N,n_feature);
    
    % Get number of rows for new_gen matrix
    % Get number of columns for new_gen matrix
    rr = size(new_gen,1);
    cc = size(new_gen,1);
    feature_choose = 1:n_feature;
    row_idx = 0; % Count up to N 
    for i = 1:rr
        row_sel = new_gen(i,:); % One row, all colums
        row_idx = row_idx + 1; 
        
        % Top row, choose col 1 and 2
        for k = feature_choose
            tempo_mat(row_idx,k) = row_sel(1,k);
        end
        
        row_idx = row_idx + 1; 
        % Botttom row, choose col 3 and 4
        for j = feature_choose
            tempo_mat(row_idx,j) = row_sel(1,(j+n_feature));
        end
    end
    mmx_reshape = tempo_mat; % Load up N by n_feature conversion
end
%---------------------------------------------------------------
