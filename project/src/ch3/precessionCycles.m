function noPrecessionCycles = precessionCycles(B0, tau, species)
% 
% Function that calculates the number of precession cycles that take place
% in the laboratory frame
% See Problem 3.3c) from "Brown, MRI"
% 
% INPUT:
%   B0 = static magnetic field strength (T)
%   tau = time of rf pulse on (ms)
%   species = 'p' proton or 'e' electron
% 
% OUTPUT:
%   noPrecessionCycles = number of precession cycles
% 
% Author: Irina Grigorescu, irina.grigorescu.15@ucl.ac.uk

% % % Prerequisites:
if nargin < 3
    species = 'p';
end

% % % % Parameters:
switch species
    case 'p'
        disp('Case a) chosen -> Proton');
        % gamma = Gyromagnetic ratio (rad/s/T)
        gamma = 2.68*1E08; % gyromagnetic ratio of proton
    case 'e'
        disp('Case b) chosen -> Electron');
        % gamma = Gyromagnetic ratio (rad/s/T)
        gamma = 658 * 2.68*1E08; % gyromagnetic ratio of electron
end

% % % % Calculate the frequency of precession
% (number of cycles per second)
freqPrecession = gamma / (2*pi) * B0;

% % % % Calculate the number of precession cycles
noPrecessionCycles = freqPrecession * (tau / 1E+03);


