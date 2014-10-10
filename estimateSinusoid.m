function [omegaList, gainList] = estimateSinusoid(y, coarseOmega, map_IfftMat,...
    map_dIfftMat, K, numStepsFine)
% INPUTS
% (i)   y  compressive measurements
% (ii)  coarseOmega coarse bin frequencies for which we have precomputed
%           compressive maps of sinusoids & their derivatives
% (iii) map_IfftMat vector representation of coarse bins 
%           map_IfftMat = S*sinusoid(omegaCoarse) (S is the comp measurement mat)
% (iv)  map_dIfftMat vector representation of derivative of signal manifold
%           at the coarse frequencies  map_dIfftMat = S*d_sinusoid(omegaCoarse)
% (v)   K  number of sinusoids in the mixture - can be inferred from noise
%           level (not implemented)
% (vi)  numStepsFine (optional) - number of steps the fine algorithm takes
%            default - 4

if ~exist('numStepsFine','var'), numStepsFine = 4;
elseif isempty(numStepsFine), numStepsFine = 4; end

% initialization
y_r = y; % residue
omegaList = [];
gainList = [];

for count = 1:K
    
    % coarse stage
    omega_add = detectNew(y_r, coarseOmega, map_IfftMat);
    
    omegaList = [omegaList; omega_add];
    gainList  = [gainList; 0]; % just a place holder 
    
    % fine stage
    [omegaList, gainList, y_r] = refineExisting(y, omegaList,...
        coarseOmega, map_IfftMat, map_dIfftMat, numStepsFine);
end