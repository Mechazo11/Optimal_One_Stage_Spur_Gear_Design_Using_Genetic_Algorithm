function new_gen = create_new_gen_spur_gears(mmx,bit_count,upper_limits,lower_limits,Nvar)
    % FUNCTION HARDCODED TO WORK WITH 3 VARIABLES ONLY    
    % upper-limit - row vector containing upper limits for Np, P and b
    % lower-limit - row vector containing upper limits for Np, P and b

    % Initialize working variables
    row_idx = size(mmx,1);
    col_idx = size(mmx,2);
    cc_mat = zeros(row_idx, col_idx); % Temporary matrix to hold children
    ccx1 = [];
    ccx2 = [];
    
    % Create (N/2) random numbers
    cross_rand = random_generator(row_idx, 0, 0.98); % (N/2) number of cross over points
    
    % For each pair of candidates, perform 'genetic mutation'
    for i = 1:row_idx
        %fprintf("Cross over iteration %d\n", i);
        in_mmx = mmx(i,:); 
        in_cross_rand = cross_rand(:,i);
        % Most important function that performs genetic cross over
        [ccx1,cxx2] = genetic_crossover(in_mmx, in_cross_rand,bit_count,upper_limits,lower_limits,Nvar);
        cons_mat = [ccx1, cxx2]; % Convert two row vector in 1
        % Create new_gen matrox
        for k = 1: col_idx
            cc_mat(i,k) = cons_mat(1,k);
        end
        new_gen = cc_mat;
    end
   
    
end

% --------------------------------------------------------------------------------------------------------------
function [child1, child2] = genetic_crossover(couple_row, rand_num, bit_count,upper_limits,lower_limits,Nvar)
    % Local function that performs the primary operation of genetic
    % crossover
    
    % Unpack upper and lower limits for each variable
    U_Np = upper_limits(1); %36
    U_P = upper_limits(2); %20
    U_b = upper_limits(3); %2.80
    
    L_Np = lower_limits(1); %12
    L_P = lower_limits(2); % 2
    L_b = lower_limits(3); % 0.450
    
    % Length of one 'chromosome' in bits
    nn_chromo = bit_count * Nvar; % if 8 bits and 2 genes then one candidate -- 16 bits 
    
    J = 2^bit_count - 1; % 65535
    
    % Probability of cross over fixed in this implementation
    crossover_prob = 0.90; % Fixed
    
    % Define ranges for bit level manipulation
    bit_range = create_bit_range_mat(bit_count, Nvar); % A matrix i.e 8 bits 2 design variables
    % [[1,2,3,4,5,6,7,8]; [9,10,11,12,13,14,15,16]]
    
    % --------------- HARDCODED SECTION FOR THREE VARIABLE ------------- %
    
    % Encode scheme
    % variable -> base_10 -> bit_count binary -> nn_chromo bit parent
    
    % Unpack row vector into six sets of gene components
    % 1 -- Np, 2 -- P, 3 -- b
    mx1 = couple_row(:,1);
    mx2 = couple_row(:,2);
    mx3 = couple_row(:,3);
    
    fx1 = couple_row(:,4);
    fx2 = couple_row(:,5);
    fx3 = couple_row(:,6);
    
    % Convert variable to X10
    mx1 = VartoX10(round(mx1), U_Np, L_Np, J); % Np 
    mx2 = VartoX10(round(mx2),U_P,L_P,J); % P
    mx3 = VartoX10(mx3,U_b,L_b,J); % b
    
    fx1 = VartoX10(round(fx1),U_Np, L_Np,J); % Ngg
    fx2 = VartoX10(round(fx2),U_P,L_P,J); % P
    fx3 = VartoX10(fx3,U_b,L_b,J); % b
    
    % Cross-over point calculation needs to be updated
    % Generate one cross over point
    N= nn_chromo - 1; 
    CrossoverIndex = cast(rand_num * N,'uint16');
    %fprintf("Cross-over point %d\n\n", CrossoverIndex);
    
    % Encode decimal to binary
    mx1 = de2bi(mx1, bit_count);
    mx2 = de2bi(mx2, bit_count);
    
    % Watchdog to prevent finite real non-negative integer error from
    % stopping the code
    try
        mx3 = de2bi(mx3, bit_count);
    catch
        warning('Finite real nonnegative integers error detected. Assigning a value of 0.5');
        mx3 = 0.5;
        mx3 = de2bi(mx3, bit_count);
    end
    
    
    fx1 = de2bi(fx1, bit_count);
    fx2 = de2bi(fx2, bit_count);
    
    try
        fx3 = de2bi(fx3, bit_count);
    catch
        warning('Finite real nonnegative integers error detected. Assigning a value of 0.5');
        fx3 = 0.5;
        fx3 = de2bi(fx3, bit_count);
    end
    
    % Make nn_chromo bit genome sequence (each design variable placed side by side)
    male_gen = [mx1, mx2, mx3]; % 3*8 bits = 24 bits
    female_gen = [fx1, fx2, fx3]; % 3*8 bits = 24 bits
    
    % Do we do crossover?
    cross_rand = random_generator(50, 0, 0.9995);
    
    % Pick one probability randomly
    rng shuffle; % Reset random seed generator
    cross_rand = randsample(cross_rand, 1, 'true', cross_rand);
    
    % Is cross rand less than
    if (crossover_prob > cross_rand)
        % Crossover, binary values
        child1 = [male_gen(1:CrossoverIndex) female_gen(CrossoverIndex+1:end)];
        child2 = [female_gen(1:CrossoverIndex) male_gen(CrossoverIndex+1:end)];
    else
        % No cross over
        child1 = male_gen;
        child2 = female_gen;
    end
    
    
    %---------------------------- Mutation operator code -----------------
    
    % Bit-Flip Mutation
    % Define ranges
    min_pm = 0.001; % Usual minimum mutation probability
    %max_pm = 1/nn_chromo; % Ref - 3 (Need to find this out)
    max_pm = 0.006; % Ref - 3
    
    % Generate nn_chromo mutation probabilities between min_pm to max_pm
    prob_mutation = random_generator(nn_chromo,min_pm,max_pm);
    % Shuffle to increase chance of randomization
    prob_mutation = randomize_array(prob_mutation);
    
    % Pick one probability randomly
    rng shuffle; % Reset random seed generator
    prob_mutation = randsample(prob_mutation, 1, 'true', prob_mutation);
    
    % Create a random probability for each bit location
    bit_mutation_prob = random_generator(nn_chromo,min_pm,max_pm);
    
    % Shuffle to increase chance of randomization
    bit_mutation_prob = randomize_array(bit_mutation_prob);
    
    % Choose a random bit location
    rng shuffle;
    bit_loc = linspace(1,nn_chromo,nn_chromo);
    bit_loc = randsample(bit_loc,1);
    bit_prob = bit_mutation_prob(1,bit_loc);
    
    % HARDCODED for two variable problem
    % Is bit_prob greater than prob_mutation? Yes, flip that bit
    if (bit_prob > prob_mutation)
        child1(1,bit_loc) = ~child1(1,bit_loc);
        child2(1,bit_loc) = ~child2(1,bit_loc);
    else
        %fprintf("No mutation!\n");
    end 
    
    %---------------------------- Mutation operator code -----------------
    % Decode scheme
    % variable <- base_10 <- bit_count binary <- nn_chromo bit parent
    % MATLAB's bit get only works on an integer
    [cx1, cx2, cx3, dx1, dx2, dx3] = X32toBin(child1, child2, bit_range);
    
%     cx1 = BintoVar(cx1, U_Np,L_Np,J);
%     cx2 = BintoVar(cx2, U_P,L_P,J);
    
    cx1 = BintoVar(cx1, U_Np,L_Np,J);
    cx2 = BintoVar(cx2, U_P,L_P,J);
    cx3 = BintoVar(cx3, U_b,L_b,J);
    
    dx1 = BintoVar(dx1, U_Np,L_Np,J);
    dx2 = BintoVar(dx2, U_P,L_P,J);
    dx3 = BintoVar(dx3, U_b,L_b,J);
    
    % Apply tool width constraint on child 1
    cx3 = apply_tooth_constraint(cx3, round(cx2));
    dx3 = apply_tooth_constraint(dx3, round(dx2));
    
    % --------------- HARDCODED SECTION FOR THREE VARIABLE ------------- %
    % Create 'children' row vectors
    child1 = [round(cx1), round(cx2), cx3]; 
    child2 = [round(dx1), round(dx2), dx3];
end
% --------------------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------------------
function b_corrected = apply_tooth_constraint(b_estimated, P_estimated)
% Applies the tooth constraint as discussed in Pg 445
b_low_lim = 9/P_estimated;
b_high_lim = 14/P_estimated;

if (b_estimated < b_low_lim & b_estimated > b_high_lim)
    %b_corrected = (b_low_lim + b_high_lim) / 2; % Take average of the two ranges
    b_corrected = b_low_lim + ((b_high_lim + b_low_lim) / 16); % Emperical equation, works
else
    b_corrected = b_estimated;
end

end

% --------------------------------------------------------------------------------------------------------------

% Dependency for genetic_crossover

% ------------------------------------------------
function [out_vector] = randomize_array(in_vector)
% Function to randomize location of elements in a row vector
% Test if in_vector is a row vector
V = isrow(in_vector);
    if (V == 0)
        in_vector = transpose(in_vector); % Make sure the in vector is a row vector
    end
rand_pos = randperm(length(in_vector)); %array of random positions
% new array with original data randomly distributed 
    for k = 1:length(in_vector)
        out_vector(k) = in_vector(rand_pos(k));
    end
% If our priginal vector was a column vector, we need to put it back in
% correct shape
    if (V == 0)
        out_vector = transpose(out_vector); % Revert it back
    end
end
% ------------------------------------------------

% ------------------------------------------------
% Xreal -- value in scale defined for the design variable
function [Xreal] = BintoVar(Xbin, U,L,J)
    Xreal = bi2de(Xbin); 
    % For some reason this has to be explicitly casted to double
    Xreal = double(Xreal);
    % Scale base-10 value to design variable scale
    Xreal = X10toVar(Xreal,U,L,J); 
end
% ------------------------------------------------

% ------------------------------------------------
function [Xreal] = X10toVar(X10,UU,LL,JJ)
  Xreal = LL + ((UU - LL)/JJ) * X10;
  % MATLAB cannot convert a float into binary
  Xreal = double(Xreal);
end
% ------------------------------------------------

% -----------------------------------------------
function [X10] = VartoX10(XVar,U,L,J)
% X10 must be fininte, non-negative integer
  
  X10 = ((XVar - L) * J)/(U - L);
  %X10 = int16(X10);
  X10 = abs(int16(X10)); % Experimental fix
end
%-----------------------------------------------

%-----------------------------------------------
function [cx1, cx2, cx3, dx1, dx2, dx3] = X32toBin(child1, child2,bit_range)
    % HARDCODED, only works for 3 variables
    % Function to automate MATLAB's way of finding and extracting binary bits
    % from an integer number
    % bit_range --> [[1,2,3,4,5,6,7,8]; [9,10,11,12,13,14,15,16], [17,18,19,20,21,22,23,24]]
    
    child1 = bi2de(child1); % bi2de converts a binary row vector b to a decimal integer.
    child2 = bi2de(child2);
    first_8 = bit_range(1,:); % First row of bit_range matrix
    second_8 = bit_range(2,:); % 2nd row, all columns
    third_8 = bit_range(3,:);
    cx1 = bitget(child1, first_8);
    cx2 = bitget(child1, second_8);
    cx3 = bitget(child1, third_8);
    
    dx1 = bitget(child2, first_8);
    dx2 = bitget(child2, second_8);
    dx3 = bitget(child2, third_8);
    
end
% -----------------------------------------------

% ----------------------------------------------
function print_genome(P)
    fprintf('Genome : [');
    fprintf('%g ', P);
    fprintf(']\n');
    P = bi2de(P);
    %fprintf("Genome in decimal --> %d\n\n", P);
end
% -----------------------------------------------

% -----------------------------------------------
function [bit_range_matrix] = create_bit_range_mat(bit_count, n_feature)
    % Function which creates a matrix containing integer locations for
    % bitwise manipulation
    nn_chromo = bit_count * n_feature;
    mmx = linspace(1,nn_chromo,nn_chromo);
    ppx = zeros(n_feature, bit_count); % Temporary matrix
    ctn = 1; % Variable to loop through each element in a row matrix
    % Use the nested loop trick to reshape
    for pp_row = 1: n_feature % Row selector
        for pp_col = 1:bit_count % column selector
            ppx(pp_row, pp_col) = mmx(1,ctn); 
            ctn = ctn + 1;
        end
    end
   bit_range_matrix = ppx; % Output
end
% -----------------------------------------------
