function [omegaList, gainList, y_r] = estimateSinusoid(y, sampledManifold, K, ...
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
%            default - 0


if ~exist('numStepsFine','var'), numStepsFine = 4;
elseif isempty(numStepsFine), numStepsFine = 4; end

N = sampledManifold.length;

if ~exist('min_sep','var'), min_sep = 0;
elseif isempty(min_sep), min_sep = 0; end

% initialization
y_r = y; % residue
omegaList = [];
gainList = [];

max_iter = 3*K; % stop adding and deleting frequencies when number of 
                % iterations exceeds this maximum

numStepsFine_intermediate = numStepsFine;
                
for iter = 1:max_iter
    
    % coarse stage
    [omega_add, gain_add, y_r] = detectNew(y_r, sampledManifold);
    
    omegaList = [omegaList; omega_add];
    gainList  = [gainList; gain_add]; % just a place holder 
    
    % refine frequencies detected so far
    [omegaList, gainList, y_r] = refineExisting(y, omegaList,...
        sampledManifold, numStepsFine_intermediate);
    
    if min_sep > 0
        % check whether two frequencies have come too close
        [omegaList, change] = pruneExisting(omegaList, min_sep);
        while change
            % refine existing frequencies to account for change
            [omegaList, gainList, y_r] = refineExisting(y, omegaList,...
                sampledManifold, numStepsFine_intermediate);
            % check again whether two frequencies have come too close
            [omegaList, change] = pruneExisting(omegaList, min_sep);
        end
    end
    
    % when we detect K frequencies we stop
    if length(omegaList) == K
        break;
    end
end

% fine stage
[omegaList, gainList, y_r] = refineExisting(y, omegaList,...
    sampledManifold, numStepsFine);
