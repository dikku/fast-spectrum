%% Define Scenario
N = 256; % Length of Sinusoid
M = 32; % Number of compressive measurements
type = 'cmplx_bernoulli'; % type of measurement matrix
                          % set to type = 'full' and M = N for
                          % non-compressive
K = 3; % # sinusoids in the mixture of sinusoids
SNR = [33 30 30]; % SNR per sinusoid (if all N measurements were made)
                  % compressive SNR scaled down roughly by a factor of M/N
                  % if isometries are met
sigma = 1;
%% Simulation Setup
NumSims = 1e4;

antidx = 1:N; % notion of starting phase

% definition of sinusoid
sinusoid    = @(omega) exp(1j*antidx(:)*omega)/sqrt(N);
% d_sinusoid  = @(omega) 1j*antidx(:)...
%     .*exp(1j*antidx(:)*omega)/sqrt(N);

S = generateMeasMat(N,M,type);
%% Algorithm parameters
overSamplingRate = 3; % Detection stage
numStepsFine     = 3; % Refinement phase
%% Algorithm preprocessing
[map_IfftMat, map_dIfftMat, coarseOmega] =...
    preProcessMeasMat(S, antidx, overSamplingRate);
%% SIMS 
omega_true = 2*pi*rand(NumSims,K);
omega_est  = zeros(NumSims,K);

gains_true = sign(randn(NumSims,K) + 1j*randn(NumSims,K))*diag(10.^(SNR/20));
gains_est  = zeros(NumSims,K);

CRB_omega = zeros(NumSims,K);

parfor sim_count = 1:NumSims
    this_gains_true = gains_true(sim_count,:);
    this_omega_true = omega_true(sim_count,:);
    
    % add measurement noise
    y_full = sigma*(randn(N,1) + 1j*randn(N,1))/sqrt(2);
    % add signal
    for count = 1:K 
        y_full = y_full +...
            this_gains_true(count)*sinusoid(this_omega_true(count));
    end
    
    % take compressive measurements
    y = S*y_full;
    
    [omegaList, gainList] = estimateSinusoid(y, coarseOmega, map_IfftMat,...
        map_dIfftMat, K, numStepsFine);
    
    omega_est(sim_count,:) = omegaList(:).';
    gains_est(sim_count,:) = gainList(:).';
    
end

%% 
omega_est_reordered  = zeros(NumSims,K);
CRB_omega_reorder    = zeros(NumSims,K);
omega_errors = zeros(NumSims,K);

wrap_2pi = @(x) angle(exp(1j*x));
parfor count = 1:NumSims
    this_omega_true = omega_true(count,:);
    this_gains_true = gains_true(count,:);
    this_omega_est = omega_est(count,:);
    
    % estimate CRB indirectly from all N measurements
    [CRB_mat, ~] = CRBAllN(this_gains_true, this_omega_true, antidx, sigma);
    this_CRB = diag(CRB_mat)*N/M; % scale CRB for compressive measurements
    this_CRB = this_CRB((2*K+1):(3*K));
    CRB_omega(count,:) = this_CRB; % we want only freq.
    
    all_possible_errors = repmat(this_omega_true,[K 1]) - ...
        repmat(this_omega_est.',[1 K]);
    [this_errors, min_idx] = min(abs(wrap_2pi(all_possible_errors)),[],1);
    
    % lazy code: this mapping may not always be a permutation, but such
    % cases will only prop up when two frequencies are very close anyway
    omega_est_reordered(count,:) = this_omega_est(min_idx); 
    omega_errors(count,:) = this_errors;
    CRB_omega_reorder(count,:) = this_CRB(min_idx);
end

%%

[f_errors,x_errors] = ecdf(omega_errors(:).^2);
[f_crb,x_crb] = ecdf(CRB_omega_reorder(:));

figure;  semilogy(10*log10(x_errors),1-f_errors,'k');
hold on; semilogy(10*log10(x_crb),1-f_crb,'r');

xlabel('Squared error in dB(\omega)','fontsize',14);
ylabel('CCDF','fontsize',14);

legend({'Algorithm','Cramer Rao Bound'});

