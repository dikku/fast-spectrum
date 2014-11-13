function [omega, gain, y_r] = detectNew(y_r, sampledManifold)

% adds a new frequency by correlating the residue with the coarse 
% templates in sampledManifold.map_IfftMat

if sampledManifold.is_eye
    R = length(sampledManifold.coarseOmega);
    N = length(y_r);
    
    possible_gains  = fft(y_r, R)/sqrt(N); 
    % phases are wrong at this point 
    % we only correct the phases 
    % for the detected frequency
else
    possible_gains = (sampledManifold.map_IfftMat'*y_r)...
        ./sampledManifold.map_IfftMat_norm_sq.';
end

[~, which_bin]  = max(abs(possible_gains));

omega = sampledManifold.coarseOmega(which_bin);
gain  = possible_gains(which_bin);

% correct phases to reflect definition of manifold
if and(sampledManifold.is_eye, sampledManifold.ant_idx(1)~=0)
    gain = gain.*exp(-1j*omega*sampledManifold.ant_idx(1));
end

% compute new residue
if sampledManifold.is_eye
    y_r = y_r - gain*exp(1j*omega*sampledManifold.ant_idx(:))/sqrt(N);
else
    y_r = y_r - gain*sampledManifold.map_IfftMat(:,which_bin);
end