function [omegaList, gainList, y_r] = estimateSinusoid(y, sampledManifold, K, ...
    numStepsFine)
% INPUTS
% (i)    y  compressive measurements
% (ii a) sampledManifold.coarseOmega coarse bin frequencies for which 
%          we have precomputed compressive maps of sinusoids & their derivatives
% (ii b) sampledManifold.map_IfftMat vector representation of coarse bins 
%           map_IfftMat = S*sinusoid(omegaCoarse) (S is the comp measurement mat)
% (ii c) sampledManifold.map_IfftMat_norm_sq: norm square of 
%           map_IfftMat = S*sinusoid(omegaCoarse) at each coarse bin
% (ii d) sampledManifold.map_dIfftMat vector representation of derivative 
%           of signal manifold at the coarse frequencies  
%           map_dIfftMat = S*d_sinusoid(omegaCoarse)
% (ii e) sampledManifold.map_d2IfftMat vector representation of the second
%           derivative of signal manifold at the coarse frequencies  
%           map_d2IfftMat = S*d2_sinusoid(omegaCoarse)
% (iii)   K  number of sinusoids to look for in the mixture
% (iv)   numStepsFine (optional) - number of steps the fine algorithm takes
%            default - 4

if ~exist('numStepsFine','var'), numStepsFine = 4;
elseif isempty(numStepsFine), numStepsFine = 4; end

% initialization
y_r = y; % residue
omegaList = [];
gainList = [];
                
for iter = 1:K
    
    % coarse stage
    [omega_add, gain_add, y_r] = detectNew(y_r, sampledManifold);
    
    omegaList = [omegaList; omega_add];
    gainList = [gainList; gain_add]; % just a place holder
    
    % recalibrate all gains jointly and update residue y_r
    [omegaList, gainList, y_r] = refineExisting(y, omegaList,...
                sampledManifold, 0);
    
    % when we detect K frequencies we stop
    if length(omegaList) == K
        break;
    end
end

% fine stage
[omegaList, gainList, y_r] = refineExisting(y, omegaList,...
    sampledManifold, numStepsFine);
