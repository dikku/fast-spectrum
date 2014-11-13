function [omegas, gains, y_r] = refineExisting(y, ...
    omegas, sampledManifold, NumFine, maxjump)
% Refines detected frequencies in omega over the [0,2pi] continuum by
% using only the coarsely sampled manifold S*x(omega_coarse) and
% S*dx(omega_coarse)/domega stored in sampledManifold
%  This is done via linearization of the manifold
% maximum jump per step is given by maxjump 
%                           (default = DFT spacing / 4 = pi/N/2)
% Complexity O(M K^2)

wrap_2pi = @(x) angle(exp(1j*x));
N = sampledManifold.length;
DFT = (2*pi/N);

if ~exist('maxjump','var'), maxjump = DFT/4;
elseif isempty(maxjump), maxjump = DFT/4; end

K = length(omegas);


if sampledManifold.is_eye
    ant_idx = sampledManifold.ant_idx(:);
    Ant_idx = sparse(1:N, 1:N, ant_idx);
    % Ant_idx2 = sparse(1:N, 1:N, ant_idx.^2);
    x_theta  = exp(1j*ant_idx*omegas.')/sqrt(N);
    dx_theta = 1j*Ant_idx...
        *exp(1j*ant_idx*omegas.')/sqrt(N);
    % d2x_theta = -Ant_idx2...
    %     *exp(1j*ant_idx*omegas.')/sqrt(N);
    
else
    delta_bin = 2*pi/length(sampledManifold.coarseOmega);
    [~,BIN_IDX] = min(abs((repmat(sampledManifold.coarseOmega(:),[1 K]) -...
        repmat(omegas(:),[1 length(sampledManifold.coarseOmega)])')),[],1);
    omegaCoarse = sampledManifold.coarseOmega(BIN_IDX)';
    omegaDelta = wrap_2pi(omegas - omegaCoarse);
    % find the nearest coarse estimate
    coarse_S_IFFT = sampledManifold.map_IfftMat(:,BIN_IDX);
    coarse_DER_S_IFFT = sampledManifold.map_dIfftMat(:,BIN_IDX);
    % d2x_theta = sampledManifold.map_d2IfftMat(:,BIN_IDX);
    % linearized version of the sinusoidal manifold around
    % the coarse estimate
    % tangent plane
    x_theta = coarse_S_IFFT + coarse_DER_S_IFFT*diag(omegaDelta);
    dx_theta = coarse_DER_S_IFFT;%  + d2x_theta*diag(omegaDelta);
end

gains = (x_theta'*x_theta)\(x_theta'*y);
y_r = y - x_theta*gains;

for iii=1:NumFine
    DER_S_IFFT_times_gain = dx_theta*diag(gains); 
    % DER2_S_IFFT_times_gain = d2x_theta*diag(gains); 
    
    % look for unexplained component along tangent in the residue
    
    % % Newton step
    % omega_jump = real(DER_S_IFFT_times_gain'*DER_S_IFFT_times_gain -...
    %     diag(DER2_S_IFFT_times_gain'*y_r))\real(DER_S_IFFT_times_gain'*y_r);
    
    % linearization about coarse estimates
    omega_jump = real(DER_S_IFFT_times_gain'*DER_S_IFFT_times_gain)\...
        real(DER_S_IFFT_times_gain'*y_r);
    
    % use the unexplained part to adjust omegas
    omegas = omegas + max(min(omega_jump, maxjump),-maxjump);
    
    if sampledManifold.is_eye
        x_theta  = exp(1j*ant_idx*omegas.')/sqrt(N);
        dx_theta = 1j*Ant_idx...
            *exp(1j*ant_idx*omegas.')/sqrt(N);
        % d2x_theta = -Ant_idx2...
        %     *exp(1j*ant_idx*omegas.')/sqrt(N);
    else
        omegaDelta = wrap_2pi(omegas - omegaCoarse);
        change_pivots = any(abs(omegaDelta) > (delta_bin/2));
        % change coarse estimate pivots if we cross boundaries when refining
        if change_pivots
            [~,BIN_IDX] = min(abs((repmat(sampledManifold.coarseOmega(:),[1 K]) -...
                repmat(omegas(:),[1 length(sampledManifold.coarseOmega)])')),[],1);
            omegaCoarse = sampledManifold.coarseOmega(BIN_IDX)';
            omegaDelta = wrap_2pi(omegas - omegaCoarse);
            
            coarse_S_IFFT = sampledManifold.map_IfftMat(:,BIN_IDX);
            coarse_DER_S_IFFT = sampledManifold.map_dIfftMat(:,BIN_IDX);
            % d2x_theta = sampledManifold.map_d2IfftMat(:,BIN_IDX);
        end
        % new template, corresponding gains and residues
        x_theta = coarse_S_IFFT + coarse_DER_S_IFFT*diag(omegaDelta);
        dx_theta = coarse_DER_S_IFFT;%  + d2x_theta*diag(omegaDelta);
    end
    
    gains = (x_theta'*x_theta)\(x_theta'*y);
    y_r = y - x_theta*gains;
end
