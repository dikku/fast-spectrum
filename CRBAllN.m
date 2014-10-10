function [CRB_mat, FIM_mat] = CRBAllN(gainList, omegaList, antidx, sigma)
% OUTPUT CRB matrix ordered as (gains(1:K), phases(1:K), frequencies(1:K))

K = length(gainList);
N = length(antidx);
sinusoid    = @(omega) exp(1j*antidx(:)*omega)/sqrt(N);
d_sinusoid  = @(omega) 1j*antidx(:)...
    .*exp(1j*antidx(:)*omega)/sqrt(N);

FIM_plane_matrix = [];

% gain portion of the CRB plane
for count = 1:K
    FIM_plane_matrix = [FIM_plane_matrix ...
        sign(gainList(count))*sinusoid(omegaList(count))];
end

% phase portion of the CRB plane
for count = 1:K
    FIM_plane_matrix = [FIM_plane_matrix ...
        1j*gainList(count)*sinusoid(omegaList(count))];
end

% frequency portion of the CRB plane
for count = 1:K
    FIM_plane_matrix = [FIM_plane_matrix ...
        gainList(count)*d_sinusoid(omegaList(count))];
end

FIM_mat = 2*real(FIM_plane_matrix'*FIM_plane_matrix)/sigma^2;
CRB_mat = inv(FIM_mat);