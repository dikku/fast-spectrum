function [omegaFine, gainEst, y_r] = polishExisting(y, omegaFine, S,...
    numsteps, maxjump)

% COSTLY REFINEMENT 
% use x(omega), dx(omega)/domega evaluated at omegaFine 
% to get close to CRB without oversampling heavily for 
% the detection phase
% Complexity O(N M K)

if ~exist('maxjump','var'), maxjump = Inf;
elseif isempty(maxjump), maxjump = Inf; end

if ~exist('numsteps','var'), numsteps = 3;
elseif isempty(numsteps), numsteps = 3; end

N = size(S,2);
M = size(S,1);
K = length(omegaFine);

% definition of sinusoid
sinusoid    = @(omega) exp(1j*(0:(N-1))'*omega)/sqrt(N);
d_sinusoid  = @(omega) 1j*(0:(N-1))'...
    .*exp(1j*(0:(N-1))'*omega)/sqrt(N);
d2_sinusoid  = @(omega) -(((0:(N-1))').^2)...
    .*exp(1j*(0:(N-1))'*omega)/sqrt(N);


omegaFine = omegaFine(:);

S_IFFT      = zeros(M,K);
DER_S_IFFT  = zeros(M,K);
DER2_S_IFFT = zeros(M,K);

for count = 1:K
    S_IFFT(:,count)      = S * sinusoid(omegaFine(count));
    DER_S_IFFT(:,count)  = S * d_sinusoid(omegaFine(count));
    DER2_S_IFFT(:,count) = S * d2_sinusoid(omegaFine(count));
end

for step_count = 1:numsteps
    
    gainEst = (S_IFFT'*S_IFFT)\(S_IFFT'*y);
    y_r = y - S_IFFT*gainEst;
    
    DER_S_IFFT_times_gain  = DER_S_IFFT*diag(gainEst);
    DER2_S_IFFT_times_gain = DER2_S_IFFT*diag(gainEst);
    
    % omega_jump = real(DER_S_IFFT_times_gain'*DER_S_IFFT_times_gain)\...
    %     real(DER_S_IFFT_times_gain'*y_r);
    
    % Newton step
    omega_jump = real(DER_S_IFFT_times_gain'*DER_S_IFFT_times_gain -...
        diag(DER2_S_IFFT_times_gain'*y_r))\real(DER_S_IFFT_times_gain'*y_r);
    
    omegaDelta = max(min(omega_jump, maxjump),-maxjump);
    omegaFine = omegaFine + omegaDelta;
    
    for count = 1:K
        S_IFFT(:,count)     = S * sinusoid(omegaFine(count));
        if step_count ~= numsteps
            DER_S_IFFT(:,count)  = S * d_sinusoid(omegaFine(count));
            DER2_S_IFFT(:,count) = S * d2_sinusoid(omegaFine(count));
        end
    end
    
end

gainEst = (S_IFFT'*S_IFFT)\(S_IFFT'*y);
y_r = y - S_IFFT*gainEst;