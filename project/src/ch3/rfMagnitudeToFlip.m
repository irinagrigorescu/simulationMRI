function B1 = rfMagnitudeToFlip(tau, theta, species)
% 
% Function that calculates the B1 magnitude of an rf pulse needed to flip
% a proton/electron spin with a certain angle
% See Problem 3.3a)b) from "Brown, MRI"
% 
% INPUT:
%   tau = time of rf pulse on (ms)
%   theta = angle to flip (deg)
%   species = 'p' proton or 'e' electron
% 
% OUTPUT:
%   B1 = the B1 field's magnitude in microT 
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


% % % % Calculate the B1 field strength in microTesla
B1 = deg2rad(theta) / (gamma * (tau/1E+03) ) * 1E+06;


