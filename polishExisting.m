function [omegas, gains, y_r] = polishExisting(y, omegas, S,...
    Numsteps, maxjump)

% COSTLY REFINEMENT 
% use x(omega), dx(omega)/domega evaluated at omegaFine 
% to get close to CRB without oversampling heavily for 
% the detection phase
% Complexity O(N M K)

N = size(S,2);
M = size(S,1);
K = length(omegas);

DFT = 2*pi/N;

if ~exist('maxjump','var'), maxjump = DFT/4;
elseif isempty(maxjump), maxjump = DFT/4; end

if ~exist('numsteps','var'), Numsteps = 3;
elseif isempty(Numsteps), Numsteps = 3; end




ant_idx = (0:(N-1)).' - (N-1)/2;

% definition of sinusoid
sinusoid    = @(omega) S * exp(1j*ant_idx*omega)/sqrt(N);
d_sinusoid  = @(omega) 1j * S * ant_idx...
    .*exp(1j*ant_idx*omega)/sqrt(N);
% d2_sinusoid  = @(omega) - S * (ant_idx.^2)...
    % .*exp(1j*ant_idx*omega)/sqrt(N);


omegas = omegas(:);

x_theta      = zeros(M,K);
dx_theta  = zeros(M,K);
% d2x_theta = zeros(M,K);

for count = 1:K
    x_theta(:,count)   = sinusoid(omegas(count));
    dx_theta(:,count)  = d_sinusoid(omegas(count));
    % d2x_theta(:,count) = d2_sinusoid(omegas(count));
end

for step_count = 1:Numsteps
    
    gains = (x_theta'*x_theta)\(x_theta'*y);
    y_r = y - x_theta*gains;
    
    x_theta_times_gain  = dx_theta*diag(gains);
    % d2x_theta_times_gain = d2x_theta*diag(gains);
    
    omega_jump = real(x_theta_times_gain'*x_theta_times_gain)\...
        real(x_theta_times_gain'*y_r);
    
    % Newton step
    % omega_jump = real(x_theta_times_gain'*x_theta_times_gain -...
    %     diag(d2x_theta_times_gain'*y_r))\real(x_theta_times_gain'*y_r);
    
    omegaDelta = max(min(omega_jump, maxjump),-maxjump);
    omegas = omegas + omegaDelta;
    
    for count = 1:K
        x_theta(:,count)       = sinusoid(omegas(count));
        if step_count ~= Numsteps
            dx_theta(:,count)  = d_sinusoid(omegas(count));
            % d2x_theta(:,count) = d2_sinusoid(omegas(count));
        end
    end
    
end

gains = (x_theta'*x_theta)\(x_theta'*y);
y_r = y - x_theta*gains;
