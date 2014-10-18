An efficient sequential algorithm to accurately estimate *continuous* frequencies/tones present in a *mixture* of sinusoids

Complexity 
==========

N - length of sinusoid
K - number of tones in the mixture
M = O(K log N) - number of compressive measurements

(i) Estimation / Cold start O(N M K)
------------------------------------

(ii) Tracking / Warm start O(M K^2) 
-----------------------------------


Overview
========

The algorithm is implemented by the "estimateSinusoid" routine which calls "detectNew" to detect new frequencies (from a coarse grid) and "refineExisting" to refine existing frequencies whenever a new frequency is detected (refinement is done on the continuum of frequencies). "pruneExisting" removes frequencies which come too close to one another and replaces them with their average (merge operation) 

The zero-padded IFFT (unit norm version) of the rows of the measurement matrix and the IFFT of the derivative manifold of the IDFT of the rows of the measurement matrix are passed as inputs to "estimateSinusoid". "preProcessMeasMat" implements the zero-padded IFFT (albeit in an inefficient manner). Can call MATLAB's IFFT routine if one desires to speed up the operation (IFFT means antidx = 0:(N-1))

"main" puts things together and benchmarks against CRB (note that when SNR is low CRB is an optimistic estimate of achievable error) ("CRBAllN" computes the CRB)


USAGE
=====

Detection
---------

N = length of each sinusoid
K = number of sinusoids in the mixture
M = number of compressive measurements
S = M by N compressive measurement matrix 
	(can be eye(N) for normal measurements)
MIN_SEP = minimum separation between sinusoids as the algorithm progresses (e.g., 1.5*2*pi/N). 

We merge two sinusoids closer than this into a single sinusoid. Should be smaller than we expect in the "truth" 

y = S*[mixture of sinusoids + noise] 
% y is compressive measurements

% precompute IFFTs upto desired oversampling rate (default is 3)
sampleManifold = preProcessMeasMat(S, N); 

% arrive at estimates from the measurements y and the sampled 
% version of the manifold in sampleManifold

[omegaEst, gainEst] = estimateSinusoid(y, sampleManifold, K, [], MIN_SEP);
