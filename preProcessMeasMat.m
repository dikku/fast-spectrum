function sampledManifold = preProcessMeasMat(S, N, overSamplingRate)

if ~exist('overSamplingRate','var'), overSamplingRate = 3;
elseif isempty(overSamplingRate), overSamplingRate = 3; end

sampledManifold.length = N;
R = round(overSamplingRate*N);
ant_idx = 0:(N-1);

%% INEFFICIENT CODE LEFT FOR READABILITY (IF ant_idx is different from 0:(N-1))

% sinusoid    = @(omega) exp(1j*ant_idx(:)*omega);
% d_sinusoid  = @(omega) 1j*ant_idx(:)...
%     .*exp(1j*ant_idx(:)*omega);
% % d2_sinusoid  = @(omega) -(freq_order(:).^2)...
%     % .*exp(1j*antidx(:)*omega);
% 
% coarseOmega = 2*pi*(0:(R-1))/(R);
% IfftMat   = zeros(N,R);
% dIfftMat  = zeros(N,R);
% % d2IfftMat = zeros(N,R);
% for count = 1:R
%     IfftMat(:,count)   = sinusoid(coarseOmega(count));
%     dIfftMat(:,count)  = d_sinusoid(coarseOmega(count));
%     % d2IfftMat(:,count) = d2_sinusoid(coarseOmega(count));
% end
% 
% 
% sampledManifold.coarseOmega = coarseOmega;
% sampledManifold.map_IfftMat =...
%     S*IfftMat; % S times x(omegaCoarse)
% sampledManifold.map_IfftMat_norm_sq = ...
%     sum(abs(sampledManifold.map_IfftMat).^2,1); 
%     % norm square of S times x(omegaCoarse)
% sampledManifold.map_dIfftMat =...
%     S*dIfftMat; % S times dx(omegaCoarse)/d omega (1st der)
% 
% % needed if we use block-Newton updates
% % sampledManifold.map_d2IfftMat =...
% %     S*d2IfftMat; 
%       % S times d^2x(omegaCoarse)/d omega^2 (2nd der)




%% EFFICIENT CODE
% Fourier transform of sensing matrix for
%  the coarse stage (can use IFFT to speed up these steps)

% SCALING OF IFFTs assumes that sinusoid(omega) takes the following form 
% sinusoid    = @(omega) exp(1j*(0:(N-1)).'*omega);

 sampledManifold.coarseOmega = 2*pi*(0:(R-1))/R; 
                               % ordering of frequencies in IFFT definition
 sampledManifold.map_IfftMat = sqrt(R*N)*ifft(S,R,2); 
                               % S times x(omegaCoarse)

 sampledManifold.map_IfftMat_norm_sq = ...
     sum(abs(sampledManifold.map_IfftMat).^2,1);
      % norm square of S times x(omegaCoarse)
 sampledManifold.map_dIfftMat = sqrt(R*N)*...
     ifft(S*sparse(1:N,1:N,1j*ant_idx),R,2); 
     % S times dx(omegaCoarse)/d omega (1st der)
     
%  sampledManifold.map_d2IfftMat = sqrt(R*N)*...
%      ifft(S*sparse(1:N,1:N,-ant_idx.^2),R,2);

