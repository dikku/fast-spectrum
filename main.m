%% Distance between phases and freq
wrap_2pi = @(x) angle(exp(1j*x));
%% Define Scenario
N = 256; % Length of Sinusoid
M = 48; % Number of compressive measurements
type = 'cmplx_bernoulli'; % type of measurement matrix
                          % set to type = 'full' and M = N for
                          % non-compressive

min_delta_omega = 4*(2*pi/N); % minimum separation between sinusoids
% SNR = 12; % SNR per sinusoid with compressive measurements
SNR = [15 12 9];% SNR per sinusoid with compressive measurements

K = length(SNR); % # sinusoids in the mixture of sinusoids
SNR_all_N = SNR + 10*log10(N/M); % SNR per sinusoid 
                                 % (if all N measurements were made)
                  
sigma = 10^(-min(SNR_all_N/20));
%% Simulation Setup
NumSims = 1e4;

antidx = 1:N; % notion of starting phase

% definition of sinusoid
sinusoid    = @(omega) exp(1j*antidx(:)*omega);
% d_sinusoid  = @(omega) 1j*antidx(:)...
%     .*exp(1j*antidx(:)*omega);

S = generateMeasMat(N,M,type);
%% Algorithm parameters
overSamplingRate = 2; % Detection stage
numStepsFine     = 3; % Refinement phase
K_est = K; % #sinusoids to look for
min_delta_omega_est = (2*pi/N); % minimum separation between 
                                  % two frequencies when we refine
%% Algorithm preprocessing
sampledManifold = preProcessMeasMat(S, antidx, overSamplingRate);
%% SIMS 
omega_true = zeros(NumSims,K);
omega_est  = zeros(NumSims,K_est);

gains_true = sigma*sign(randn(NumSims,K) + 1j*randn(NumSims,K))*...
    diag(10.^(SNR_all_N/20));
gains_est  =  zeros(NumSims,K_est);

parfor sim_count = 1:NumSims
    
    this_gains_true = gains_true(sim_count,:);
    
    this_omega_true = 2*pi*rand(1,1);
    while length(this_omega_true) < K
        temp = 2*pi*rand(1,1);
        if min(abs(wrap_2pi(this_omega_true-temp))) > min_delta_omega
            this_omega_true = [this_omega_true temp];
        end
    end
    omega_true(sim_count,:) = this_omega_true;
    
    % add measurement noise
    y_full = sigma*(randn(N,1) + 1j*randn(N,1))/sqrt(2);
    % add signal
    for count = 1:K 
        y_full = y_full +...
            this_gains_true(count)*sinusoid(this_omega_true(count));
    end
    
    % take compressive measurements
    y = S*y_full;
    
    [omegaList, gainList] = estimateSinusoid(y, sampledManifold, K_est,...
        numStepsFine, min_delta_omega_est);
    
    omegaListFull = Inf*ones(K_est,1);
    omegaListFull(1:length(omegaList)) = omegaList(:);
    
    gainListFull = Inf*ones(K_est,1);
    gainListFull(1:length(gainList)) = gainList(:);
    
    omega_est(sim_count,:) = omegaListFull.';
    gains_est(sim_count,:) = gainListFull.';
    
end

%% 
omega_est_reordered  = zeros(NumSims,K);
omega_errors = zeros(NumSims,K);
CRB_omega = zeros(NumSims,K);

parfor count = 1:NumSims
    this_omega_true = omega_true(count,:);
    this_gains_true = gains_true(count,:);
    this_omega_est = omega_est(count,:);
    
    % estimate CRB indirectly from all N measurements
    [CRB_mat, ~] = CRBAllN(this_gains_true, this_omega_true, antidx, sigma);
    this_CRB = diag(CRB_mat)*N/M; % scale CRB for compressive measurements
    this_CRB = this_CRB((2*K+1):(3*K));
    CRB_omega(count,:) = this_CRB; % we want only freq.
    
    all_possible_errors = repmat(this_omega_true,[K_est 1]) - ...
        repmat(this_omega_est.',[1 K]);
    [this_errors, min_idx] = min(abs(wrap_2pi(all_possible_errors)),[],1);
    
    % This mapping may not always be a permutation, but such
    % cases will come up when two frequencies are very close
    omega_est_reordered(count,:) = this_omega_est(min_idx); 
    omega_errors(count,:) = this_errors;
end

%% Plot results

[f_errors,x_errors] = ecdf(omega_errors(:).^2);
[f_crb,x_crb] = ecdf(CRB_omega(:));

f1 = figure;  semilogy(10*log10(x_errors),1-f_errors,'k');
hold on; semilogy(10*log10(x_crb),1-f_crb,'r');

xlabel('Squared frequency estimation error in dB','fontsize',16);
ylabel('CCDF','fontsize',16);

legend({'Algorithm','Cramer Rao Bound'},'fontsize',16,'location','southwest');
% print('-depsc2','-r300','./Results.eps');

%% Plot one - by - one
f2 = figure;
for count = 1:NumSims
    this_omega_true = omega_true(count,:);
    this_omega_est = omega_est(count,:);
    this_errors = omega_errors(count,:);
    this_gains_true = gains_true(count,:);
    this_gains_est = gains_est(count,:);
    this_CRB = CRB_omega(count,:);
    
    if (max(abs(this_errors).^2./this_CRB)) > 20
        figure(f2);
        hold off;
        stem(this_omega_true, abs(this_gains_true),'r-o');
        hold on;
        plot(this_omega_est, abs(this_gains_est),'kx');
        pause(1);
    end
end
close(f2);



