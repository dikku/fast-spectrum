function [omega_add, gain_add, y_r] = detectNew(y, sampledManifold)

% adds a new frequency by correlating the residue with the coarse 
% templates in sampledManifold.map_IfftMat

xcorr = abs(y'*sampledManifold.map_IfftMat)...
    ./sampledManifold.map_IfftMat_norm_sq;
[~, which_bin]  = max(xcorr);

omega_add = sampledManifold.coarseOmega(which_bin);
gain_add  = (sampledManifold.map_IfftMat(:,which_bin)' * y)/...
    sampledManifold.map_IfftMat_norm_sq(which_bin);

y_r = y - gain_add*sampledManifold.map_IfftMat(:,which_bin);