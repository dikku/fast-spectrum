A gradient based, sequential algorithm to detect frequencies/tones in a mixture of sinusoids

Complexity O(N M K) (N - length of sinusoid, M - number of compressive measurements, K - number of tones in the mixture)

The algorithm is implemented by the "estimateSinusoid" routine which calls "detectNew" to detect new frequencies (from a coarse grid) and "refineExisting" to refine existing frequencies whenever a new frequency is detected (refinement is done on the continuum of frequencies). "estimateSinusoid" assumes that the zero-padded IFFT of the measurement matrix and the IFFT of the derivative manifold are computed upfront (this depends only on the measurement matrix). "preProcessMeasMat" implements the zero-padded IFFT (albeit in an inefficient manner). Can call MATLAB's FFT routine if one desires

"main" puts things together and benchmarks against CRB (note that when SNR is low CRB is an optimistic estimate of achievable error)


USAGE:

N = length of your sinusoid
M = number of compressive measurements
S = M by N compressive measurement matrix
y = compressive measurements

% define antenna indexing
antidx = 1:N; 

% precompute IFFT's upto desired oversampling rate (default is 3)
[map_IfftMat, map_dIfftMat, coarseOmega] = preProcessMeasMat(S, 1:N); 

% arrive at estimates from the measurements y and the information about the measurement setup stored in coarseOmega, map_IfftMat, map_dIfftMat

[omegaEst, gainEst] = estimateSinusoid(y, coarseOmega, map_IfftMat, map_dIfftMat, K, numStepsFine);