function [omegaList, gainList] = estimateSinusoid(y, sampledManifold, K, ...
    numStepsFine, min_sep)
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
% (iii)   K  number of sinusoids in the mixture - can be inferred from noise
%           level (not implemented)
% (iv)   numStepsFine (optional) - number of steps the fine algorithm takes
%            default - 4
% (v)    min_sep (optional) - we will merge two sinusoids that are closer
%              than this in frequency
%            default - pi/N


if ~exist('numStepsFine','var'), numStepsFine = 4;
elseif isempty(numStepsFine), numStepsFine = 4; end

N = sampledManifold.length;

if ~exist('min_sep','var'), min_sep = pi/N;
elseif isempty(min_sep), min_sep = pi/N; end

% initialization
y_r = y; % residue
omegaList = [];
gainList = [];

max_iter = 3*K; % stop adding and deleting frequencies when number of 
                % iterations exceeds this maximum

for iter = 1:max_iter
    
    % coarse stage
    omega_add = detectNew(y_r, sampledManifold);
    
    omegaList = [omegaList; omega_add];
    gainList  = [gainList; 0]; % just a place holder 
    
    % fine stage
    [omegaList, gainList, y_r] = refineExisting(y, omegaList,...
        sampledManifold, numStepsFine);
    
    % check whether two frequencies have come too close
    [omegaList, change] = pruneExisting(omegaList, min_sep);
    while change
        
        % if change = true -- we need to refine again as omegaList has been
        % pruned
        if change
            [omegaList, gainList, y_r] = refineExisting(y, omegaList,...
                sampledManifold, numStepsFine);
        end
        
        % check again whether two frequencies have come too close
        [omegaList, change] = pruneExisting(omegaList, min_sep);
    end
    
    % when we detect K frequencies we stop
    if length(omegaList) == K
        break;
    end
end
