function freq = resonanceFrequency(B0, gamma)
%
% Function that calculates the frequency associated 
%  with the RF field at a given static magnetic field value
% 
% Input:
%    B0    = External magnetic field (Tesla)
%    gamma = Gyromagnetic ratio (rad/s/T)
% Output:
%    freq  = Calculated frequency (MHz)
%
% Author: Irina Grigorescu, irina.grigorescu.15@ucl.ac.uk

% Parameters:
if nargin < 2
    gamma = 2.68*1E08; % gyromagnetic ratio
end

% Frequency
freq = gamma .* B0 ./ (2*pi) .* 1E-06; % transforms values to MHz

