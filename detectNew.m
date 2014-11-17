function [omega, gain, y_r] = detectNew(y_r, sampledManifold)

% adds a new frequency by correlating the residue with the coarse 
% templates in sampledManifold.map_IfftMat

if sampledManifold.is_eye
    R = length(sampledManifold.coarseOmega);
    N = length(y_r);
    
    % energy = 1; % for all frequencies
    possible_gains  = fft(y_r, R)/sqrt(N); 
    
    % phases are wrong at this point 
    
    % correct gains
    possible_gains = possible_gains.*...
        exp(- 1j*sampledManifold.ant_idx(1) *...
        sampledManifold.coarseOmega(:));
    
else
    energy = sampledManifold.map_IfftMat_norm_sq.';
    possible_gains = (sampledManifold.map_IfftMat'*y_r)./energy;
end

% choose the frequency which reduces the residue by the largest amount
if sampledManifold.is_eye
    [~, which_bin]  = max((abs(possible_gains).^2));
else
    [~, which_bin]  = max((abs(possible_gains).^2) .* energy);
end

omega = sampledManifold.coarseOmega(which_bin);
% gain  = possible_gains(which_bin);
% 
% % compute new residue
% if sampledManifold.is_eye
%     y_r = y_r - gain*exp(1j*omega*sampledManifold.ant_idx(:))/sqrt(N);
% else
%     y_r = y_r - gain*sampledManifold.map_IfftMat(:,which_bin);
% end

% refine the detected coarse frequency to prevent "leakage"
[omega, gain, y_r] = refine_new_freq(y_r, omega, sampledManifold);

end

function [omega, gain, y_r] = refine_new_freq(y, omega, sampledManifold)
wrap_2pi = @(x) angle(exp(1j*x));

if sampledManifold.is_eye
    N = length(y);
    ant_idx = sampledManifold.ant_idx(:);
    x_theta  = exp(1j*ant_idx*omega)/sqrt(N);
    dx_theta = 1j*ant_idx...
        .*exp(1j*ant_idx*omega)/sqrt(N);
    d2x_theta = -(ant_idx.^2)...
        .*exp(1j*ant_idx*omega)/sqrt(N);
    energy = 1;
else
    [~,bin] = min(abs(wrap_2pi(sampledManifold.coarseOmega - omega)));
    omegaCoarse = sampledManifold.coarseOmega(bin);
    omegaDelta = wrap_2pi(omega - omegaCoarse);
    coarse_S_IFFT = sampledManifold.map_IfftMat(:,bin);
    coarse_DER_S_IFFT = sampledManifold.map_dIfftMat(:,bin);
    d2x_theta = sampledManifold.map_d2IfftMat(:,bin);
    x_theta = coarse_S_IFFT + coarse_DER_S_IFFT*omegaDelta;
    dx_theta = coarse_DER_S_IFFT + d2x_theta*omegaDelta;
    energy = x_theta'*x_theta;
end


% % NON-COHERENT
% want derivatives of abs(x_theta'*y)^2/energy
% r = log(energy);
if sampledManifold.is_eye
    % energy = constant ( = 1 )
    % d_r = 0;
    % d2_r = 0;
    
    g =  x_theta'*y;
    d_g = dx_theta'*y;
    d2_g = d2x_theta'*y;
else
    d_r = 2*real(x_theta'*dx_theta)/energy;
    d2_r = 2*(real(x_theta'*d2x_theta) + ...
        (dx_theta'*dx_theta))/energy - d_r^2;
    
    g =  x_theta'*y/sqrt(energy);
    d_g = (dx_theta'*y)/sqrt(energy) - g*d_r/2;
    d2_g = (d2x_theta'*y)/sqrt(energy) - ...
        (d_g + g*d_r/4)*d_r - g*d2_r/2;
end

% f = abs(g)^2;
d_f = 2*real(g'*d_g);
d2_f = 2*real(g'*d2_g) + 2*abs(d_g)^2;

der1 = - d_f;
der2 = - d2_f;


if der2 > 0
    omega_next = omega - der1/der2;
else
    DFT = 2*pi/sampledManifold.length;
    omega_next = omega - max(min(der1, DFT/4),-DFT/4);
end

omega = omega_next;

if sampledManifold.is_eye
    x_theta  = exp(1j*ant_idx*omega)/sqrt(N);
else
    [~,bin_new] = min(abs(wrap_2pi(...
        sampledManifold.coarseOmega - omega)));
    
    changed_pivots = (bin_new ~= bin);
    if changed_pivots
        bin = bin_new;
        omegaCoarse = sampledManifold.coarseOmega(bin);
        omegaDelta = wrap_2pi(omega - omegaCoarse);
        coarse_S_IFFT = sampledManifold.map_IfftMat(:,bin);
        coarse_DER_S_IFFT = sampledManifold.map_dIfftMat(:,bin);
    else
        omegaDelta = wrap_2pi(omega - omegaCoarse);
    end
    
    x_theta = coarse_S_IFFT + coarse_DER_S_IFFT*omegaDelta;
    energy = (x_theta'*x_theta);
end

gain = (x_theta'*y)/energy;
y_r = y - gain*x_theta;

end