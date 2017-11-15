function B0 = fieldStrength(freq, gamma)
%
% Function that calculates the field strength associated 
%  with the frequency of an RF field for a given spin species.
% 
% Input:
%    freq  = RF frequency (MHz)
%    gamma = Gyromagnetic ratio (rad/s/T)
% Output:
%    B0    = External magnetic field (Tesla)
%
% Author: Irina Grigorescu, irina.grigorescu.15@ucl.ac.uk

% Parameters:
if nargin < 2
    gamma = 2.68*1E08; % gyromagnetic ratio
end

% Field Strength
B0 = freq * (2*pi) / gamma * 1E06; % in Tesla