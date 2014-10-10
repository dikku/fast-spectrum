function [omegaFine, gainEst, y_r] = refineExisting(y, omegaFine,...
    CoarseOmegaBins, map_IfftMat, map_dIfftMat, numStepsFine, maxjump)
% CRUX OF THE ALGORITHM

if ~exist('maxjump','var'), maxjump = Inf;
elseif isempty(maxjump), maxjump = Inf; end

K = length(omegaFine);
delta_bin = 2*pi/length(CoarseOmegaBins);

% find the nearest coarse estimate
[~,BIN_IDX] = min(abs((repmat(CoarseOmegaBins(:),[1 K]) -...
    repmat(omegaFine(:),[1 length(CoarseOmegaBins)])')),[],1);

omegaCoarse = CoarseOmegaBins(BIN_IDX)';
omegaDelta = omegaFine - omegaCoarse;
coarse_S_IFFT = map_IfftMat(:,BIN_IDX);
coarse_DER_S_IFFT = map_dIfftMat(:,BIN_IDX);

% linearized version of the sinusoidal manifold around the coarse estimate
% tangent plane
S_template = coarse_S_IFFT + coarse_DER_S_IFFT*diag(omegaDelta);
gainEst = S_template\y;
y_r = y - S_template*gainEst; % compute the residue

for iii=1:numStepsFine
    d_omega_MAT = coarse_DER_S_IFFT*diag(gainEst); 
    % look for unexplained component along tangent in the residue
    omega_jump = (d_omega_MAT'*d_omega_MAT)\real(d_omega_MAT'*y_r);
    % use the unexplained part to adjust omegas
    omegaFine = omegaFine + min(omega_jump, maxjump);
    omegaDelta = omegaFine - omegaCoarse;
    
    % change coarse estimate pivots if we cross boundaries when refining
    if any(abs(omegaDelta) > (delta_bin/2)) 
        [~,BIN_IDX] = min(abs((repmat(CoarseOmegaBins(:),[1 K]) -...
            repmat(omegaFine(:),[1 length(CoarseOmegaBins)])')),[],1);
        
        omegaCoarse = CoarseOmegaBins(BIN_IDX)';
        omegaDelta = omegaFine - omegaCoarse;
        coarse_S_IFFT = map_IfftMat(:,BIN_IDX);
        coarse_DER_S_IFFT = map_dIfftMat(:,BIN_IDX);
    end
    
    % new template, corresponding gains and residues
    S_template = coarse_S_IFFT + coarse_DER_S_IFFT*diag(omegaDelta);
    gainEst = S_template\y;
    y_r = y - S_template*gainEst;
end