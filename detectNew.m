function omega_add = detectNew(y_r, sampledManifold)

% adds a new frequency by correlating the residue with the coarse 
% templates in sampledManifold.map_IfftMat

xcorr = abs(y_r'*sampledManifold.map_IfftMat)...
    ./sampledManifold.map_IfftMat_norm_sq;
[~, which_bin]  = max(xcorr);

omega_add = sampledManifold.coarseOmega(which_bin);