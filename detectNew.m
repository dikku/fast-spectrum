function omega_add = detectNew(y_r, omegaCoarse, map_IfftMat)

xcorr    = abs(y_r'*map_IfftMat)./sqrt(sum(abs(map_IfftMat).^2,1));
[~, which_bin]  = max(xcorr);
omega_add = omegaCoarse(which_bin);