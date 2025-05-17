% Define parameter values for the HHM (Hidden Hodgkin-Huxley Model)
% These represent conductances, reversal potentials, time constants, etc.
keys = ["gL" "gK" "gNa" "gCa" "gT" "gAHP" "vL" "vK" "vNa" "vCa" ...
    "vGS" "tau_h1" "tau_n1" "tau_r1" "tau_h0" "tau_n0" "tau_r0" ...
    "fi_n" "fi_h" "fi_r" "k1" "kCa" "eps" ...
    "alpha" "beta" "tetm" "teth" "tetn" "tetr" "teta" "tetb" "tets" ...
    "sigm" "sigh" "sign" "sigr" "siga" "sigb" "sigs" "tetg" "tetHg" ...
    "tet_th" "tet_tn" "tet_tr" "sig_th" "sig_tn" "sig_tr" "sigHg"];
values = [2.25 45.0 37.5 0.5 0.5 9.0 -60.0 -80.0 55.0 140.0 -85.0 ...
    500.0 100.0 17.5 1.0 1.0 40.0 0.75 0.75 0.2 15.0 22.5 3.75e-5 ...
    5.0 1.0 -30.0 -39.0 -32.0 -67.0 -63.0 0.4 -39.0 ...
    15.0 -3.1 8.0 -2.0 7.8 -0.1 8.0 30.0 -39.0 ...
    -57.0 -80.0 68.0 -3.0 -26.0 -2.2 8.0];
params = dictionary(keys, values);

%% Define Network Size and Synaptic Strengths
num_neurons = 100;  % Number of neurons in each population
inhibition_strength = 1;  % (placeholder) Not yet used
excitation_strength = 1;  % (placeholder) Not yet used

% Define number of connections per neuron based on literature ratios
Ngs = round(6 * num_neurons / 244);  % GPe → STN
Ngg = round(40 * num_neurons / 244); % GPe → GPe
Nsg = 1;  % STN → GPe

%% Initialize Synaptic Connectivity Matrices
GPe_to_GPe = zeros(num_neurons);  % GPe inhibiting GPe
GPe_to_STN = zeros(num_neurons);  % GPe inhibiting STN
STN_to_GPe = zeros(num_neurons);  % STN exciting GPe

%% Build Synaptic Connections (Manual Permutation Approach)

for i = 1:num_neurons
    % Random inhibitory targets for GPe → STN and GPe → GPe
    gpe_stn = randperm(num_neurons, Ngs);  % GPe to STN
    gpe_gpe = randperm(num_neurons, Ngg);  % GPe to GPe
    
    % Random excitatory target from STN → GPe
    stn_gpe = randperm(num_neurons, Nsg);
    
    % Apply connections
    GPe_to_GPe(i, gpe_gpe) = 1;
    GPe_to_STN(i, gpe_stn) = 1;
    STN_to_GPe(i, stn_gpe) = 1;
end

%% Avoid Self-Connections in GPe → GPe
for i = 1:num_neurons
    if GPe_to_GPe(i, i) ~= 0
        GPe_to_GPe(i, i) = 0;  % Remove self-connection
        
        % Add a new random connection to maintain degree
        available = find(GPe_to_GPe(i, :) == 0);
        new_target = available(randi(length(available)));
        GPe_to_GPe(i, new_target) = 1;
    end
end

%% Set Activation Matrices for Model
% These will be used in the simulation code later

STN_GPe = STN_to_GPe;  % STN → GPe
GPe_STN = GPe_to_STN;  % GPe → STN
GPe_GPe = GPe_to_GPe;  % GPe → GPe

%% Initialize Neural State Vector
% Initial membrane potentials; -70 mV is typical resting potential
Data_zero = zeros(12 * num_neurons, 1);
for i = 1:num_neurons
    Data_zero(i) = -70;                % For STN or first half
    Data_zero(6 * num_neurons + i) = -70;  % For GPe or second half
end

%% Example Applied Values (for simulation)
Iapp1 = -5.2;  % Applied current to neurons
GS = 2.5;      % Gain or conductance scaling (model-specific)

%% 📌 Suggested Apviews (Conceptual Highlights)
% - Shows how **circuit-level structure** (STN ↔ GPe) is implemented.
% - Demonstrates the **balance of excitation/inhibition** via sparse matrices.
% - Highlights **manual control** over topology vs. random network assumptions.
% - This structure can be modified to explore **dysfunction in Parkinsonian loops**.

