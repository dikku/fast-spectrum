function sampledManifold = preProcessMeasMat(S, overSamplingRate)

if ~exist('overSamplingRate','var'), overSamplingRate = 3;
elseif isempty(overSamplingRate), overSamplingRate = 3; end

M = size(S,1);
N = size(S,2);

sampledManifold.length = N;
R = round(overSamplingRate*N);

% good choice of initial phase for linearization (needed to reduce
% overSamplingRate
ant_idx = (0:(N-1)) - (N-1)/2; 

if M == N
    sampledManifold.is_eye = norm(eye(N) - S,'fro') == 0;
else
    sampledManifold.is_eye = false;
end

%% INEFFICIENT CODE - BUT AN INTERPRETATION OF THE IFFT OPERATION

% sinusoid    = @(omega) exp(1j*ant_idx(:)*omega)/sqrt(N);
% d_sinusoid  = @(omega) 1j*ant_idx(:)...
%     .*exp(1j*ant_idx(:)*omega)/sqrt(N);
% % d2_sinusoid  = @(omega) -(ant_idx(:).^2)...
%     % .*exp(1j*ant_idx(:)*omega)/sqrt(N);
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
% % needed when we use block-Newton updates
% % sampledManifold.map_d2IfftMat =...
% %     S*d2IfftMat; 
% %     % S times d^2x(omegaCoarse)/d omega^2 (2nd der)




%% EFFICIENT CODE - USE IFFTs
% Fourier transform of sensing matrix for
%  the coarse stage (can use IFFT to speed up these steps)

% SCALING OF IFFTs assumes that sinusoid(omega) takes the following form 
% sinusoid    = @(omega) exp(1j*(0:(N-1)).'*omega);

sampledManifold.coarseOmega = 2*pi*(0:(R-1))/R;  % omegaCoarse
                         % ordering of frequencies in IFFT's definition
sampledManifold.ant_idx = ant_idx;

if ~sampledManifold.is_eye
    
    sampledManifold.map_IfftMat = R/sqrt(N)*ifft(S,R,2)*...
        diag(exp(1j*sampledManifold.coarseOmega*ant_idx(1)));
    % S times x(omegaCoarse)
    
    sampledManifold.map_IfftMat_norm_sq = ...
        sum(abs(sampledManifold.map_IfftMat).^2,1);
    % norm square of S times x(omegaCoarse)
    sampledManifold.map_dIfftMat = R/sqrt(N)*...
        ifft(S*sparse(1:N,1:N,1j*ant_idx),R,2)*...
        diag(exp(1j*sampledManifold.coarseOmega*ant_idx(1)));
    % S times dx(omegaCoarse)/d omega (1st der)
    
    sampledManifold.map_d2IfftMat = R/sqrt(N)*...
        ifft(S*sparse(1:N,1:N,-ant_idx.^2),R,2)*...
        diag(exp(1j*sampledManifold.coarseOmega*ant_idx(1)));
    % S times d^2x(omegaCoarse)/d omega^2 (2nd der)

end
end