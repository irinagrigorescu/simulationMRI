function spexc = excessSpins(B0, T, gamma)
%
% Function that calculates the spin excess as a fraction of N 
%  for protons at a given external magnetic field value (tesla)
% 
% Input:
%    B0    = External magnetic field (Tesla)
%    T     = Temperature (Kelvin)
%    gamma = Gyromagnetic ratio(rad/s/T)
% Output:
%    spexc = Calculated spin excess (as a fraction)
%
% Author: Irina Grigorescu, irina.grigorescu.15@ucl.ac.uk
%                           irinagry@gmail.com

% Parameters:
reducedPlanck  = 1.05*(10^-34); % Reduced Planck constant (J*s)
k              = 1.38*(10^-23); % Boltzmann's constant (J/K)

if nargin < 2
    T              = 300;       % Temperature (K)
end

if nargin < 3
    gamma          = 2.68*1E08; % Gyromagnetic ratio
end

% Angular frequency
omega0 = gamma*B0;

% Spin excess
spexc = (reducedPlanck*omega0) / (2*k*T);