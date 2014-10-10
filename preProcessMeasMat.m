function [map_IfftMat, map_dIfftMat, coarseOmega] =...
    preProcessMeasMat(S, antidx, overSamplingRate)
if ~exist('overSamplingRate','var'), overSamplingRate = 3;
elseif isempty(overSamplingRate), overSamplingRate = 3; end

N = length(antidx);
sinusoid    = @(omega) exp(1j*antidx(:)*omega)/sqrt(N);
d_sinusoid  = @(omega) 1j*antidx(:)...
    .*exp(1j*antidx(:)*omega)/sqrt(N);


coarseOmega = 2*pi*(0:(N*overSamplingRate-1))/(N*overSamplingRate);
IfftMat  = zeros(N,N*overSamplingRate);
dIfftMat = zeros(N,N*overSamplingRate);
for count = 1:(N*overSamplingRate)
    IfftMat(:,count)     = sinusoid(coarseOmega(count));
    dIfftMat(:,count)     = d_sinusoid(coarseOmega(count));
end

% Fourier transform of sensing matrix for the coarse stage 
% (naively implemented here)
map_IfftMat  = S*IfftMat;
map_dIfftMat = S*dIfftMat;