An efficient sequential algorithm to accurately estimate *continuous* frequencies/tones present in a *mixture* of sinusoids


Complexity 
==========

N - length of sinusoid
K - number of tones in the mixture
M = O(K log N) - number of compressive measurements

(i) Estimation / Cold start O(N M K)
-------------------------------------

(ii) Tracking / Warm start O(M K^2) 
------------------------------------


Overview
========

The algorithm is implemented by the "estimateSinusoid" routine which calls "detectNew" to detect new frequencies (from a coarse grid) and "refineExisting" to refine existing frequencies whenever a new frequency is detected (refinement is done on the continuum of frequencies). "pruneExisting" removes frequencies which come too close to one another and replaces them with their average (merge operation) 

The zero-padded IFFT (unit norm version) of the rows of the measurement matrix and the IFFT of the derivative manifold of the IDFT of the rows of the measurement matrix are passed as inputs to "estimateSinusoid". "preProcessMeasMat" implements the zero-padded IFFT

"main" puts things together and benchmarks against CRB (note that when SNR is low CRB is an optimistic estimate of achievable error) ("CRBAllN" computes the CRB)


USAGE
=====

N = length of each sinusoid
K = number of sinusoids in the mixture
M = number of compressive measurements
S = M by N compressive measurement matrix 
	(can be eye(N) for normal measurements)
% y is compressive measurements
y = S*[mixture of sinusoids + noise] or S*[mixture of sinusoids] + noise

% precompute IFFTs upto desired 
% oversampling rate (default is 3) 
sampledManifold = preProcessMeasMat(S, N); 


(i) Estimation / Cold start O(N M K)
------------------------------------

MIN_SEP = minimum separation between sinusoids 
		  as the algorithm progresses (e.g., 1.5*2*pi/N). 

% We merge two sinusoids closer than MIN_SEP in 
% frequency into a single sinusoid
% MIN_SEP Should be smaller than the frequency separation 
% among frequencies in the "truth" 



% arrive at estimates from the measurements y and the sampled 
% version of the manifold in sampleManifold

[omegaEst, gainEst] = estimateSinusoid(y, sampledManifold, K, [], MIN_SEP);

(ii) Tracking / Warm start O(M K^2) 
-----------------------------------

omegaPastEstimate = An approximate estimate of the frequencies of the K sinusoids in y; for example from the past time-step (tracking)

[omegaRefined, gainRefined, y_residue] = refineExisting(y, omegaPastEstimate, sampledManifold);