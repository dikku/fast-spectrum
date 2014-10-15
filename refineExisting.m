function [omegaFine, gainEst, y_r] = refineExisting(y, omegaFine,...
    sampledManifold, numStepsFine, maxjump)
% Refines detected frequencies in omegaFine over the [0,2pi] continuum by
% using only the coarsely sampled manifold S*x(omega_coarse) and
% S*dx(omega_coarse)/domega stored in sampledManifold
%  This is done via linearization of the manifold

if ~exist('maxjump','var'), maxjump = Inf;
elseif isempty(maxjump), maxjump = Inf; end

K = length(omegaFine);
delta_bin = 2*pi/length(sampledManifold.coarseOmega);

% find the nearest coarse estimate
[~,BIN_IDX] = min(abs((repmat(sampledManifold.coarseOmega(:),[1 K]) -...
    repmat(omegaFine(:),[1 length(sampledManifold.coarseOmega)])')),[],1);
omegaCoarse = sampledManifold.coarseOmega(BIN_IDX)';
omegaDelta = omegaFine - omegaCoarse;
coarse_S_IFFT = sampledManifold.map_IfftMat(:,BIN_IDX);
coarse_DER_S_IFFT = sampledManifold.map_dIfftMat(:,BIN_IDX);
% coarse_DER2_S_IFFT = sampledManifold.map_d2IfftMat(:,BIN_IDX);

% linearized version of the sinusoidal manifold around the coarse estimate
% tangent plane
S_template = coarse_S_IFFT + coarse_DER_S_IFFT*diag(omegaDelta);
gainEst = S_template\y;
y_r = y - S_template*gainEst; % compute the residue

for iii=1:numStepsFine
    DER_S_IFFT_times_gain = coarse_DER_S_IFFT*diag(gainEst); 
    % DER2_S_IFFT_times_gain = coarse_DER2_S_IFFT*diag(gainEst); 
    
    % look for unexplained component along tangent in the residue
    
    % Newton step
    % omega_jump = real(DER_S_IFFT_times_gain'*DER_S_IFFT_times_gain -...
    %     diag(DER2_S_IFFT_times_gain'*y_r))\real(DER_S_IFFT_times_gain'*y_r);
    
    % linearization about coarse estimates
    omega_jump = real(DER_S_IFFT_times_gain'*DER_S_IFFT_times_gain)\...
        real(DER_S_IFFT_times_gain'*y_r);
    
    % use the unexplained part to adjust omegas
    omegaFine = omegaFine + min(omega_jump, maxjump);
    omegaDelta = omegaFine - omegaCoarse;
    
    % change coarse estimate pivots if we cross boundaries when refining
    if any(abs(omegaDelta) > (delta_bin/2))
        [~,BIN_IDX] = min(abs((repmat(sampledManifold.coarseOmega(:),[1 K]) -...
            repmat(omegaFine(:),[1 length(sampledManifold.coarseOmega)])')),[],1);
        omegaCoarse = sampledManifold.coarseOmega(BIN_IDX)';
        omegaDelta = omegaFine - omegaCoarse;
        coarse_S_IFFT = sampledManifold.map_IfftMat(:,BIN_IDX);
        coarse_DER_S_IFFT = sampledManifold.map_dIfftMat(:,BIN_IDX);
        % coarse_DER2_S_IFFT = sampledManifold.map_d2IfftMat(:,BIN_IDX);
    end
    
    % new template, corresponding gains and residues
    S_template = coarse_S_IFFT + coarse_DER_S_IFFT*diag(omegaDelta);
    gainEst = S_template\y;
    y_r = y - S_template*gainEst;
end