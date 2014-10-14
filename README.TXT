A gradient based, sequential algorithm to detect frequencies/tones in a mixture of sinusoids

Complexity O(N M K) (N - length of sinusoid, M - number of compressive measurements, K - number of tones in the mixture)

The algorithm is implemented by the "estimateSinusoid" routine which calls "detectNew" to detect new frequencies (from a coarse grid) and "refineExisting" to refine existing frequencies whenever a new frequency is detected (refinement is done on the continuum of frequencies). "estimateSinusoid" assumes that the zero-padded IFFT of the measurement matrix and the IFFT of the derivative manifold are computed upfront (this depends only on the measurement matrix). "preProcessMeasMat" implements the zero-padded IFFT (albeit in an inefficient manner). Can call MATLAB's FFT routine if one desires

"main" puts things together and benchmarks against CRB (note that when SNR is low CRB is an optimistic estimate of achievable error)


USAGE:

N = length of each sinusoid
K = number of sinusoids in the mixture
M = number of compressive measurements
S = M by N compressive measurement matrix 
	% (can be eye(N) for normal measurements)
min_sep = minimum separation between sinusoids in the mixture ( e.g., 2*pi/N)
y = compressive measurements

% define antenna indexing
antidx = 1:N; 

% precompute IFFTs upto desired oversampling rate (default is 3)
sampleManifold = preProcessMeasMat(S, antidx); 

% arrive at estimates from the measurements y and the sampled 
% version of the manifold in sampleManifold

[omegaEst, gainEst] = estimateSinusoid(y, sampleManifold, K, [], min_sep);