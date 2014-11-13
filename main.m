%% Distance between phases and freq
wrap_2pi = @(x) angle(exp(1j*x));
%% Define Scenario
N = 256; % Length of Sinusoid
DFT = (2*pi/N); % DFT spacing
M = N; % Number of compressive measurements
type = 'full'; % type of measurement matrix
               % set to type = 'full' and M = N for
               % non-compressive

% minimum separation between sinusoids
min_delta_omega = 3*DFT; 
max_delta_omega = Inf;

% Account for phase wrap around effects
min_delta_omega = min(min_delta_omega, pi);
max_delta_omega = min(max_delta_omega, 2*pi - 2*min_delta_omega);

% effective SNR per measurement with compressive 
% measurements for the K sinusoids in the mixture
% SNR = 12; 
K = 5;
SNR = 33 * ones(1,K);

K = length(SNR); % # sinusoids in the mixture of sinusoids
SNR_all_N = SNR + 10*log10(N/M); % actual SNR per measurement
                  
sigma = 10^(-min(SNR_all_N/20));
%% Simulation Setup
NumSims = 1e4;

% definition of sinusoid
ant_idx = (0:(N-1)).';
sinusoid    = @(omega) exp(1j*ant_idx*omega)/sqrt(N);
% d_sinusoid  = @(omega) 1j*ant_idx...
%     .*exp(1j*(0:(N-1))'*omega)/sqrt(N);

S = generateMeasMat(N,M,type);
%% Algorithm parameters
overSamplingRate = 2; % Detection stage
numStepsFine     = 6; % Refinement phase
K_est = K; % #sinusoids to look for
%% Algorithm preprocessing
sampledManifold = preProcessMeasMat(S, overSamplingRate);
%% SIMS 
omega_true = zeros(NumSims,K);
omega_est  = zeros(NumSims,K_est);
residue = zeros(NumSims,1);

gains_true = sigma*sign(randn(NumSims,K) + 1j*randn(NumSims,K))*...
    diag(10.^(SNR_all_N/20));
gains_est  =  zeros(NumSims,K_est);

parfor sim_count = 1:NumSims
    
    this_gains_true = gains_true(sim_count,:);
    
    this_omega_true = 2*pi*(rand(1,1)-0.5);
    
    while length(this_omega_true) < K
        temp = this_omega_true(end) + min_delta_omega +...
            (max_delta_omega-min_delta_omega)*rand(1);
        this_omega_true = [this_omega_true wrap_2pi(temp)];
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
    
    [omegaList, gainList, y_r] = estimateSinusoid(y, sampledManifold,...
        K_est, numStepsFine);
    
    % [omegaList, gainList, y_r] = polishExisting(y, omegaList, S,...
    %     numStepsFine);
    
    residue(sim_count) = y_r'*y_r;
    
    omegaListFull = Inf*ones(K_est,1);
    omegaListFull(1:length(omegaList)) = omegaList(:);
    
    gainListFull = Inf*ones(K_est,1);
    gainListFull(1:length(gainList)) = gainList(:);
    
    omega_est(sim_count,:) = omegaListFull.';
    gains_est(sim_count,:) = gainListFull.';
    
end
% account for antenna indexing changes
gains_est = gains_est.*...
        exp(1j*(sampledManifold.ant_idx(1) - ant_idx(1))*omega_est);
    
%% COMPUTE ESTIMATION ERRORS AND BENCHMARK AGAINST CRAMER RAO BOUND

omega_est_reordered  = zeros(NumSims,K);
omega_errors = zeros(NumSims,K);
CRB_omega = zeros(NumSims,K);

parfor count = 1:NumSims
    
    this_omega_true = omega_true(count,:);
    this_gains_true = gains_true(count,:);
    this_omega_est = omega_est(count,:);
    
    % estimate CRB indirectly from all N measurements
    [CRB_mat, ~] = CRBAllN(this_gains_true, this_omega_true, N, sigma);
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

[f_errors,x_errors] = ecdf(abs(omega_errors(:)).^2/DFT^2);
[f_crb,x_crb] = ecdf(CRB_omega(:)/DFT^2);

f1 = figure;  semilogy(10*log10(x_errors),1-f_errors,'k');
hold on; semilogy(10*log10(x_crb),1-f_crb,'r');
semilogy(20*log10([1 1]),[1/NumSims 1],'b');

% xlabel('Squared frequency estimation error in dB (relative to DFT))',...
%     'fontsize',16);
ylabel('CCDF','fontsize',16);

legend({'Squared Error','Cramer Rao Bound','DFT spacing'},...
    'fontsize',16,'location','southwest');

[f_res,x_res] = ecdf(residue/(M-K_est)/sigma^2);
% [f_res,x_res] = ecdf(residue/M/sigma^2);
f2 = figure; semilogy(x_res, 1-f_res);
% xlabel('Residue relative to noise level');
ylabel('CCDF');

%% DEBUG PLOT: Displays large error scenarios case by case

f3 = figure;
count_err = 0;
for count = 1:NumSims
    this_omega_true = omega_true(count,:);
    this_omega_est = omega_est(count,:);
    this_errors = omega_errors(count,:);
    this_gains_true = gains_true(count,:);
    this_gains_est = gains_est(count,:);
    this_CRB = CRB_omega(count,:);
    
    if (max(abs(this_errors).^2)/max(this_CRB)) > 100
        % (max(abs(this_errors).^2./this_CRB)) > 100
        figure(f3);
        hold off;
        stem(this_omega_true, abs(this_gains_true),'r-o');
        hold on;
        plot(this_omega_est, abs(this_gains_est),'kx');
        % pause(1);
        count_err = count_err + 1;
    end
end
disp(count_err);
close(f3);