function [omegaList, change] = pruneExisting(omegaList, min_sep)
% removes frequencies which are closer than min_sep and replaces them with
% their mean

wrap_2pi = @(x) angle(exp(1j*x));
pdiffMat = @(x) repmat(x(:),[1 length(x)]) - repmat(x(:),[1 length(x)]).';

pdistMat = abs(wrap_2pi(pdiffMat(omegaList))) +...
        diag(Inf*ones(size(omegaList)));
very_close = pdistMat < min_sep;

change = any(very_close(:));

if any(very_close(:))
    
    K = length(omegaList);
    
    [I,J,dist_IJ] = find(pdistMat.*very_close);
    
    [~,IDX] = min(dist_IJ);
    
    i = I(IDX);
    j = J(IDX);
    
    mu = (omegaList(i) + omegaList(j))/2;
    omegaList(i) = mu;
    
    omegaList = omegaList([1:(j-1) (j+1):K]);
    
end