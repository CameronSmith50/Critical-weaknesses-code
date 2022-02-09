% This code is designed to simulate the SEIRD model with high and low risk
% populations both within the community and withing care homes in order to
% simulate the effect of sheilding.
%
% Any commented numbers to the right of values are the default values.
%
% Original author: Kit Yates
% Edited by: Cameron Smith
%
% First created: 24/11/20
% Last modified: 31/01/22

clear
close all

% Start the timer
tic

% Initialise the save name
save_name='herd_immunity';

% Decide the size of the population 
N = 1e6; % Default: 1e6

% Decide how many peope are infected initially
I_init = 10;  % Default: 10

% Decide the proportion of higher risk individuals
h=0.07;  % Default: 0.07

% Decide the fraction of those that live in the main community
c=0.897; % Default: 0.897

% Decide how many long term care facilities we would like. We split these
% into small, medium and large
n_small = 120; % Default: 120
n_med = 48; % Default: 48
n_large = 24; % Default: 24
n_vec = [n_small, n_med, n_large];
n = n_small + n_med + n_large;

% Size of each care-home
s_small = 20;  % Default: 20
s_med = 50;  % Default: 50
s_large = 100;  % Default: 100
s_vec = [s_small, s_med, s_large];

% Number of people in each group, and update N. First is the low-risk,
% second is the high-risk in the community and the remaining are the
% numbers in care homes, ordered, small then medium then large
N_vec = round([N*(1-h),N*h*c,s_small*ones(1,n_small),s_med*ones(1,n_med),s_large*ones(1,n_large)]);
N = sum(N_vec);
h = sum(N_vec(2:end))/N;
c = N_vec(2)/sum(N_vec(2:end));

% Calculate how many subpopulations in total
T = n+2;

% Define the transition rate from E to I
sigma = 1/5;  % Default: 1/5

% Define the transition rate from I to R or D
Gamma = 1/2;  % Default: 1/2

% Define Beta 0 the transmission probability per contact
beta_0 = 3*Gamma;  % Default: 3*Gamma

% Define the mortality risk in the low risk group
alpha_L = 0.001;  % Default: 0.001

% Define the mortality rate in the high risk group
alpha_H = 0.05;  % Default: 0.

% Rate of waning immunity (set to 8 months)
nu = 0;  % Default: 0

% Rate of external infection. Note this only affects the two community
% subpopulations
eta_vec = 0*[1,1,zeros(1,n)];  % Default: 0 for all

% Define the intra-carehome contact rate
delta = 0.9;  % Default: 0.9

% Asymptomatic proportion
prop_asymp = 1/3;  % Default: 1/3

% Define the increase factor in beta0 caused by symptomatic transmission
% compared to asymptomatic transmission, and the altered IFRs
split_inf = 0;  % Set to 1 if you wish for the infected class to be split.
symp_trans_increase = 1.1;
asymp_trans_decrease = 0.8;
if split_inf == 1
    alpha_L = 1/(1-prop_asymp)*alpha_L;
    alpha_H = 1/(1-prop_asymp)*alpha_H;
end

% Define epsilon the threshold below which sheilding is relaxed. This value
% is the number of symptomatic cases in the population per day.
epsilon = 60/100000*N;  % Default: 60/100000*N

% Create the name for saving
save_name = sprintf('%s_threshold=%03d_',save_name,round(epsilon*100000/N));
if sum(eta_vec) > 0
    save_name = [save_name,'ext_'];
end

% Give a range of different sheilding options. Select which you want by
% setting its value to 1. Only one can be selected at a time.
NS=1;       % No Shielding
IS=0;       % Imperfect shielding
PS=0;       % Perfect shielding

% Check that we haven't chosen too many Shielding types
if NS+IS+PS > 1
    error('You have chosen too many shielding types')
end
% Or too few
if NS+IS+PS < 1
    error('You have not chosen enough shielding types')
end

% No shielding
if NS == 1
    % Define the sheilding parameters
    q_s1 = 1;  % Low risk community to low risk community. Default: 1
    q_s2 = 1;  % High risk community to high risk community. Default: 1
    q_s3 = 1;  % Care home to care homes. Default: 1
    q_s4 = 1;  % Low risk to high risk in the community and vice versa. Default: 1
    q_s5 = 1;  % Low risk community to Carehome and vice versa. Default: 1
    q_s6 = 1;  % High risk community to Carehome and vicee versa. Default: 1
    
    % Add a descriptor to the save name
    save_name = [save_name,'NS'];
end

% Imperfect shielding
if IS == 1
    % Define the sheilding parameters
    q_s1 = 1;  % Low risk community to low risk community. Default: 1
    q_s2 = 0.2;  % High risk community to high risk community. Default: 0.2
    q_s3 = 0.2;  % Care home to care homes. Default: 0.2
    q_s4 = 0.2;  % Low risk to high risk in the community and vice versa. Default: 0.2
    q_s5 = 0.2;  % Low risk community to Carehome and vice versa. Default: 0.2
    q_s6 = 0;  % High risk community to Carehome and vicee versa. Default: 0
    
    % Add a descriptor to the save name
    save_name = [save_name,'IS'];
end

% Perfect shielding
if PS == 1
    % Define the sheilding parameters
    q_s1 = 1;  % Low risk community to low risk community. Default: 1
    q_s2 = 0;  % High risk community to high risk community. Default: 0
    q_s3 = 0;  % Care home to care homes. Default: 0
    q_s4 = 0;  % Low risk to high risk in the community and vice versa. Default: 0
    q_s5 = 0;  % Low risk community to Carehome and vice versa. Default: 0
    q_s6 = 0;  % High risk community to Carehome and vicee versa. Default: 0
    
    
    % Add a descriptor to the save name
    save_name = [save_name,'PS'];
end

% Initialise Q the shielding matrix
Q = zeros(T);

% Based on these values, define the matrix Q which moderates the interaction
% matrix when shielding
Q(1,1) = q_s1;  % Low risk community to low risk community
Q(2,2) = q_s2;  % High risk community to high risk community
Q(3:T,3:T) = q_s3;  % Care home to care homes
Q(2,1) = q_s4;  % Low risk to high risk in the community and vice versa
Q(1,2) = q_s4;  % Low risk to high risk in the community and vice versa
Q(3:T,1) = q_s5;  % Low risk community to Carehome and vice versa
Q(1,3:T) = q_s5;  % Low risk community to Carehome and vice versa
Q(3:T,2) = q_s6;  % High risk community to Carehome and vice versa
Q(2,3:T) = q_s6;  % High risk community to Carehome and vice versa

% Define the probability of transmission from the low risk community into a
% carehome
P_HF_L = (1-h)*(1-delta)/(1-h+c*h);

% Define the probability of transmission from the high risk community into a
% carehome
P_HF_HC = (1-delta-P_HF_L);

% Define the probability of tansmission from the carehome into the low risk
% community - remember to divide by n as there are n care homes
P_L_HF = P_HF_L/((1-h));

% Define the probability of transmission from the carehome into the high risk
% community
P_HC_HF = P_HF_HC/(c*h);

% Divide through the proportion of contacts for community to each carehome by the number of
% carehomes
P_L_HF = P_L_HF.*N_vec(3:end)/N;
P_HC_HF = P_HC_HF.*N_vec(3:end)/N;

% Define the probability of tranmission in the low risk community
P_L_L = (1-h)*(1-sum(P_L_HF))/(1-h+c*h);

% Define the probability of transmission from high risk into low risk
% community
P_L_HC = c*h*(1-sum(P_L_HF))/(1-h+c*h);

% Define the probability of transmission from low risk into high risk
% community
P_HC_L = (1-h)*(1-sum(P_HC_HF))/(1-h+c*h);

% Define the probability of transmission from the high risk to the high risk
% community
P_HC_HC = c*h*(1-sum(P_HC_HF))/(1-h+c*h);

% Define the matrix of transition probabilities
beta_mat = zeros(T,T);

% Define the matrix of within community transmission
delta_vec = delta*ones(n,1);
delta_mat = diag(delta_vec);

% Write this to beta_mat
beta_mat(3:T,3:T) = delta_mat;

% Define the transmission matrix beta
% Define the transition from low risk community to low risk community
beta_mat(1,1) = P_L_L;
% Define the transition from high risk community to high risk community
beta_mat(2,2) = P_HC_HC;
% Define the transition from low risk community to high risk community
beta_mat(2,1) = P_HC_L;
% Define the transition from high risk community to low risk community
beta_mat(1,2) = P_L_HC;
% Define the transition from Low risk community to care homes
beta_mat(3:T,1) = P_HF_L;
% Define the transition from High risk community to care homes
beta_mat(3:T,2) = P_HF_HC;
% Define the transition from care homes to low risk comunity
beta_mat(1,3:T) = P_L_HF;
% Define the transition from care homes to high risk communities
beta_mat(2,3:T) = P_HC_HF;
% Define the transition from care homes to care homes
beta_mat(3:T,3:T) = delta_mat;

% Scale beta_mat by the appropriate rate
beta_mat = beta_0*beta_mat;

% Define the vector which specifies the mortality rate
alpha_vec = [alpha_L;alpha_H*ones(n+1,1)];

% Define the final time we are interested in simulating until (in days)
T_final = 600;  % Default: 600

% Define the recording time interval
rec_step = 1;  % Default: 1

% Define the time vector to plot
time_vector = 0:rec_step:T_final;

% Calculate the number of recording steps required
num_rec_steps = T_final/rec_step;

% Define M the number of repeats we will undertake
M = 10;  % Default: 100

% Define the vector which will reocrd the mean cell numbers at each time
% point and for each repeat
S_mean = zeros(T,num_rec_steps+1,M);  % Susceptible
E_mean = zeros(T,num_rec_steps+1,M);  % Exposed
if split_inf == 1
    IS_mean=zeros(T,num_rec_steps+1,M); % Infected symptomatic
    IA_mean=zeros(T,num_rec_steps+1,M); % Infected asymptomatic
else
    I_mean = zeros(T,num_rec_steps+1,M);  % Infected
end
R_mean = zeros(T,num_rec_steps+1,M);  % Recovered
D_mean = zeros(T,num_rec_steps+1,M);  % Dead

% Define the shielding time matrix
shielding_mean = zeros(2,M);

% Define the matrices for the number of new cases
new_cases_per_day = zeros(num_rec_steps+1,M);  % Number of new cases per day
new_cases_per_day_group = zeros(T,num_rec_steps+1,M);  % Number of new cases per day by group

% Proportion initially immune
prop_immune = 0.0;
immune_L = prop_immune*N_vec(1);
immune_HC = prop_immune*N_vec(2);

% Initialise the mean numbers at the first time_step
% Susceptibles
S_mean(1,1,:) = round(N*(1-h)-I_init-immune_L);  % Non-vulnerable susceptibles minus the initial immune and initial infected
S_mean(2,1,:) = round(N*h*c)-immune_HC;  % Vulnerable susceptibles in the community minus initial immune
S_mean(3:T,1,:) = N_vec(3:T)'*ones(1,M);  % LTC member suscebtibles

% Infecteds
if split_inf == 1
    IS_mean(1,1,:) = I_init;  % Non-vulnerable initial infected symptomatic
else
    I_mean(1,1,:) = I_init;  % Non-vulnerable initial infected
end

% Recovered
R_mean(1,1,:) = immune_L;  % Non-vulnerable initially immune
R_mean(2,1,:) = immune_HC;  % Vulneranble in the commmunity initially immune

% Add the initial cases into the new_cases_per_day
new_cases_per_day(1,:) = I_init;

% Indicator for the console
fprintf('\n')
fprintf('Running code for %d repeats:\n',M)
fprintf('\n')
fprintf('Percentage  0                                                100\n')
fprintf('            |                                                 |\n')
fprintf('            |') 

%Run through the for loop for repeats
for m=1:M
    
    % Add to console if needed
    if mod(m/M*100,2) < 1e-6 || abs(mod(m/M*100,2)-2) < 1e-6
        fprintf('|')
    end
    
    %Define the initial time
    t = 0;
    
    %Initialise the times to help recording
    t_after = t;
   
    %Initialise a vector for each disease status which tells us how many
    %people in each subpopulation
    S = zeros(T,1); % Susceptible
    E = zeros(T,1); % Exposed
    if split_inf == 1
        IS = zeros(T,1);  % Infected symptomatic
        IA = zeros(T,1);  % Infected asymptomatic
    else
        I = zeros(T,1); % Infected
    end
    R = zeros(T,1); % Recovered
    D = zeros(T,1); % Dead
    new_cases=zeros(T,1); % Counter for new cases in each subpopulation
    
    % Initialise the state vectors
    % Susceptibles
    S(1) = round(N*(1-h)-I_init-immune_L);  % Non-vulnerable susceptibles minus the initial immune and initial infected
    S(2) = round(N*h*c)-immune_HC;  % Vulnerable susceptibles in the community minus initial immune
    S(3:T) = N_vec(3:T);  % LTC member suscebtibles
    
    % Infecteds
    if split_inf == 1
        IS(1) = I_init;  % Non-vulnerable initial infected symptomatic
    else
        I(1) = I_init;  % Non-vulnerable initial infected
    end
    
    % Recovered
    R(1) = immune_L;  % Non-vulnerable initially immune
    R(2) = immune_HC;  % Vulneranble in the commmunity initially immune    
    
    % Predefine the propensity function vector
    a = zeros(5*T+4*n+2,1);
    
    % Start off shielding
    SHIELDING = 1;
    shielding_lifted = 0;
    
    % Write a variable for if the number of cases is decreasing
    DECREASING = 0;
    
    % Define beta_mat_a and Q_mat_A to the the transition matrix we are 
    % actually using (i.e. modified for shielding)
    beta_mat_a = beta_mat.*Q;
    Q_mat_a = Q;
    eta_vec_a = zeros(1,length(eta_vec));
    
    % Enter the while loop to loop through time
    while t<T_final
        
        % Calculate the propensity functions
        if split_inf == 1
            a(1:T) = eta_vec_a'.*S + S.*symp_trans_increase.*IS.*diag(beta_mat_a)./N_vec' + S.*IA.*symp_trans_increase.*diag(beta_mat_a)./N_vec';  % S+I -> E+I for all groups
        else
            a(1:T) = eta_vec_a'.*S + S.*I.*diag(beta_mat_a)./N_vec';  % S+I -> E+I for all groups
        end
        a((T+1):(2*T)) = sigma*E;  % E -> I for all groups
        if split_inf == 1
            a(2*T+1) = Gamma*(1-alpha_vec(1))*IS(1) + Gamma*IA(1);  % I -> R for non-vulnerable
            a((2*T+2):(3*T)) = Gamma.*(1-alpha_vec(2:T)).*IS(2:T) + Gamma*IA(2:T);  % I -> R for all vulnerable
            a(3*T+1) = Gamma.*alpha_vec(1)*IS(1);  % I -> D for non vulnerable
            a((3*T+2):(4*T)) = Gamma.*alpha_vec(2:T).*IS(2:T);  % I -> D for all vulnerable
        else
            a(2*T+1) = Gamma*(1-alpha_vec(1))*I(1);  % I -> R for non-vulnerable
            a((2*T+2):(3*T)) = Gamma.*(1-alpha_vec(2:T)).*I(2:T);  % I -> R for all vulnerable
            a(3*T+1) = Gamma.*alpha_vec(1)*I(1);  % I -> D for non vulnerable
            a((3*T+2):(4*T)) = Gamma.*alpha_vec(2:T).*I(2:T);  % I -> D for all vulnerable
        end
        a((4*T+1):(5*T)) = nu*R;  % R -> S for all groups
        
        %Now specificy the cross infection rates
        if split_inf == 1
            a((5*T+1):(5*T+n))=beta_mat_a(3:T,1)*symp_trans_increase.*IS(1).*S(3:T)./N_vec(1) + beta_mat_a(3:T,1)*symp_trans_increase.*IA(1).*S(3:T)./N_vec(1);  % Low risk community to carehomes
            a((5*T+n+1):(5*T+2*n))=beta_mat_a(3:T,2)*symp_trans_increase.*IS(2).*S(3:T)./N_vec(2) + beta_mat_a(3:T,2)*symp_trans_increase.*IA(2).*S(3:T)./N_vec(2);  % High risk community to carehomes
            a((5*T+2*n+1):(5*T+3*n))=beta_mat_a(1,3:T)'.*symp_trans_increase.*IS(3:T)*S(1)./N_vec(3:T)' + beta_mat_a(1,3:T)'.*symp_trans_increase.*IA(3:T)*S(1)./N_vec(3:T)';  % Carehomes to low risk comunity
            a((5*T+3*n+1):(5*T+4*n))=beta_mat_a(2,3:T)'.*symp_trans_increase.*IS(3:T)*S(2)./N_vec(3:T)' + beta_mat_a(2,3:T)'.*symp_trans_increase.*IA(3:T)*S(2)./N_vec(3:T)';  % Carehomes to high risk community
            a(5*T+4*n+1)=beta_mat_a(2,1)*symp_trans_increase*IS(1)*S(2)./N_vec(1) + beta_mat_a(2,1)*symp_trans_increase.*IA(1)*S(2)./N_vec(1);  % Low risk community to high risk community
            a(5*T+4*n+2)=beta_mat_a(1,2)*symp_trans_increase*IS(2)*S(1)./N_vec(2) + beta_mat_a(1,2)*symp_trans_increase.*IA(2)*S(1)./N_vec(2);  % High risk community to low risk community
        else
            a((5*T+1):(5*T+n)) = beta_mat_a(3:T,1)*I(1).*S(3:T)./N_vec(1);  % Low risk community to carehomes
            a((5*T+n+1):(5*T+2*n)) = beta_mat_a(3:T,2)*I(2).*S(3:T)./N_vec(2);  % High risk community to carehomes
            a((5*T+2*n+1):(5*T+3*n)) = beta_mat_a(1,3:T)'.*I(3:T)*S(1)./N_vec(3:T)';  % Carehomes to low risk comunity
            a((5*T+3*n+1):(5*T+4*n)) = beta_mat_a(2,3:T)'.*I(3:T)*S(2)./N_vec(3:T)';  % Carehomes to high risk community
            a(5*T+4*n+1) = beta_mat_a(2,1)*I(1)*S(2)./N_vec(1);  % Low risk community to high risk community
            a(5*T+4*n+2) = beta_mat_a(1,2)*I(2)*S(1)./N_vec(2);  % High risk community to low risk community
        end
        
        % Calculate the sum of the propensity functions
        a0 = sum(a);
        
        % Determine the time for the next reaction
        tau = (1/a0)*log(1/rand);
        
        % Update the time
        t = t+tau;
        
        % Generate a random number which will allow us to decide which
        % reaction is going to take place next
        r1 = a0*rand;
        
        % Define cum_sum to be the first propensity initially
        cum_sum=a(1);
        
        % Calculate the checkpoint to augment using the quickfind mex file
        j = quickfind(a,1,r1);
        
        % Previous state values
        S_before = sum(S);
        E_before = sum(E);
        
        % Implement the reaction which corresponds to j
        if j <= T
            % Remove a susceptible of type j and add an exposed of type j
            S(j) = S(j)-1;
            E(j) = E(j)+1;
            
        elseif j <= 2*T
            % Find the appropriate index
            j = j-T;
            % Move an E of type j to an I of type j
            E(j) = E(j)-1;
            if split_inf == 1 && rand < prop_asymp
                IA(j) = IA(j)+1;
            elseif split_inf == 1
                IS(j) = IS(j)+1;
            else
                I(j) = I(j)+1;
            end
            new_cases(j) = new_cases(j)+1;  % We have a new case so add to counter
            
        elseif j <= 3*T
            % Find the appropriate index
            j = j-2*T;
            % Move an I of type j to an R of type j
            if split_inf == 1 && rand*(IA(j)+IS(j)*(1-alpha_vec(j))) < IA(j)
                IA(j) = IA(j)-1;
            elseif split_inf == 1
                IS(j) = IS(j)-1;
            else
                I(j) = I(j)-1;
            end
            R(j) = R(j)+1;
            
        elseif j <= 4*T
            % Find the appropriate index
            j = j-3*T;
            % Move an I of type j to a D of type j
            if split_inf == 1
                IS(j) = IS(j)-1;
            else
                I(j) = I(j)-1;
            end
            D(j) = D(j)+1;

        elseif j <= 5*T
            % Find the appropriate index
            j = j-4*T;
            % Move an R of type j to an S of type j
            R(j) = R(j)-1;
            S(j) = S(j)+1;
            
        elseif j <= 5*T+n
            % Find the appropriate index
            j = j-5*T+2;
            % Remove a susceptible of type j and add an exposed of type j
            S(j) = S(j)-1;
            E(j) = E(j)+1;
            
        elseif j <= 5*T+2*n
            % Find the appropriate index
            j = j-5*T-n+2;
            % Remove a susceptible of type j and add an exposed of type j
            S(j) = S(j)-1;
            E(j) = E(j)+1;
            
        elseif j <= 5*T+3*n
            % Remove a susceptible of type j and add an exposed of type j
            S(1) = S(1)-1;
            E(1) = E(1)+1;
            
        elseif j <= 5*T+4*n
            % Remove a susceptible of type j and add an exposed of type j
            S(2) = S(2)-1;
            E(2) = E(2)+1;
            
        elseif j == 5*T+4*n+1
            % Remove a susceptible of type j and add an exposed of type j
            S(2) = S(2)-1;
            E(2) = E(2)+1;
            
        elseif j == 5*T+4*n+2
            % Remove a susceptible of type j and add an exposed of type j
            S(1) = S(1)-1;
            E(1) = E(1)+1;
        end
        
        % Calculate the times for recording
        t_before = t_after;
        t_after = t;
        
        % Calculate the indices of the time step before and the current time
        % step in terms of recording
        ind_before = ceil((t_before+eps)/rec_step);
        ind_after = min(floor(t_after/rec_step),num_rec_steps);
        
        % Find out how many time-steps to write to
        steps_to_write = ind_after-ind_before+1;
        
        % If we have crossed one or more recording time-steps and we have not reached the final
        % time
        if steps_to_write > 0
            S_mean(:,ind_before+1:ind_after+1,m) = S_mean(:,ind_before+1:ind_after+1,m)+S*ones(1,steps_to_write); % Susceptible
            E_mean(:,ind_before+1:ind_after+1,m) = E_mean(:,ind_before+1:ind_after+1,m)+E*ones(1,steps_to_write); % Exposed
            if split_inf == 1
                IS_mean(:,ind_before+1:ind_after+1,m) = IS_mean(:,ind_before+1:ind_after+1,m)+IS*ones(1,steps_to_write); % Infected symptomatic
                IA_mean(:,ind_before+1:ind_after+1,m) = IA_mean(:,ind_before+1:ind_after+1,m)+IA*ones(1,steps_to_write); % Infected asymptomatic
            else
                I_mean(:,ind_before+1:ind_after+1,m) = I_mean(:,ind_before+1:ind_after+1,m)+I*ones(1,steps_to_write); % Infected
            end
            R_mean(:,ind_before+1:ind_after+1,m) = R_mean(:,ind_before+1:ind_after+1,m)+R*ones(1,steps_to_write); % Recovered
            D_mean(:,ind_before+1:ind_after+1,m) = D_mean(:,ind_before+1:ind_after+1,m)+D*ones(1,steps_to_write); % Dead;
            
            % Store these new cases
            new_cases_per_day(ind_before+1:ind_after+1,m) = new_cases_per_day(ind_before+1:ind_after+1,m) + sum(new_cases)/steps_to_write;
            new_cases_per_day_group(:,ind_before+1:ind_after+1,m) = new_cases_per_day_group(:,ind_before+1:ind_after+1,m) + new_cases*ones(1,steps_to_write)/steps_to_write;

            % Reset the new_cases vector
            new_cases = zeros(T,1);
            
            % Ensure we have at least seven days of data. If so, find the
            % number of cases in that week
            if ind_after > 6
                cases_per_week = sum(new_cases_per_day((ind_after-5):(ind_after+1)));
            end
            
            % Gradient of new cases
            new_cases_grad = (new_cases_per_day(ind_after+1,m) - new_cases_per_day(ind_before,m))/(steps_to_write*rec_step);
            
            % Decide if we are now decreasing by firstly checking if we
            % weren't before, and secondly checking if the gradient is
            % significantly negative. Note that the threshold for this may
            % not always work and may need tweaking depending on if you
            % expect a flatter initial curve (need a smaller number such as
            % -80). Cannot be 0 due to stochastic fluctuations so needs to 
            % be large enough in absolute value to avoid such fluctuations 
            % triggering this. 
            if DECREASING == 0 && new_cases_grad < -120
                DECREASING = 1;
            end

            % Find out if shielding should end
            if SHIELDING == 1 && ind_after > 6 && DECREASING == 1 && cases_per_week < epsilon/(1-prop_asymp)
                % Update so that we are no longer shielding
                SHIELDING = 0;
                % Record the time at which Shielding was lifted
                shielding_lifted = t;
                % Reset the transition matrix to its original value
                beta_mat_a = beta_mat;
                % Update the Q matrix
                Q_mat_a = eye(T);
                eta_vec_a = eta_vec;
                % Add the shielding time
                shielding_mean(2,m) = shielding_lifted;
            end
        end
        
    end
    
end

% Finish the timer
TIME = toc;
fprintf('     Time taken = %.4f\n\n',TIME)

% Append to the save_name
save_name = [save_name, '.mat'];

% Save the data to the relvant save_name
save(save_name);
