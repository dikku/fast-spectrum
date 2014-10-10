function S = generateMeasMat(N,M,type)
% INPUTS
% N is the length of sinusoid
% M is the number of random projections
% type (optional) type of measurement matrix
%       SUPPORTED TYPES: complex gaussian : cmplx_gauss (default)
%                        real gaussian    : gauss
%                        bernoulli {pm 1} : bernoulli
%                        complex bernoulli {pm 1,pm j} : cmplx_bernoulli
%                        full             : All N (one by one, not compressive)

if ~exist('type','var'), type = 'cmplx_gauss';
elseif isempty(type), type = 'cmplx_gauss'; end

S = (randn(M,N) + 1j*randn(M,N));

switch type
    case 'cmplx_gauss'
        S = S/sqrt(2*N);
    case 'gauss'
        S = real(S)/sqrt(N);
    case 'bernoulli'
        S = exp(1j*floor(angle(S)/pi)*pi)/sqrt(N);
    case 'cmplx_bernoulli'
        S = exp(1j*floor(angle(S)/(pi/2))*pi/2)/sqrt(N);
    case 'full'
        if M ~= N
            error('Number of measurements M must be equal to the length of the sinusoid N');
        end
        S = eye(N);
    otherwise
        error('Measurement ensemble not supported');
end
