% % % % IRINA GRIGORESCU
% % % % DATE: 28-Feb-2017
% % % % CHAPTER 18.1.4 Generating a constant transverse magnetisation
% % % % 
% % % % This script plots the flip angle as a function of
% % % % rf pulse number in order to force a constant
% % % % transverse magnetisation from one rf pulse to the next
% % % % 

clear all; close all; clc
addpath ../../helpers

%% Loading some data
run config; params = ans; clear ans;

%% The relaxation function
E12 = @(tn,t12) exp(-tn./t12);

%%  Function that calculates the angle which depends on the previous angle
% % theta_i = asin [sin(theta_0) / ...
% %           (E1 * sin(theta_0) / tan(theta_i-1) + 1 - E1) ]
Thetai = @(theta0, thetaim1, TR, T1) real( asin (sin(theta0) ./ ...
          (E12(TR, T1) .* sin(theta0)./tan(thetaim1) + (1 - E12(TR,T1)))));
      
                  
%% Plot figure
TR     =  10; %(ms)
T1     = 950; %(ms) GM
theta0 =  12;

n = 100; % 100 iterations
thetaValues    = zeros(1,n);
thetaValues(1) = theta0;

% Calculate the angle iteratively knowing theta0
for idx = 2:n
    
    % Calculate iteratively the next value for the 
    % rf flip angle knowing that theta_0 = 13 deg,
    % TR = 10ms and T1 = 950ms
    thetaValues(idx) = rad2deg(...
        Thetai(deg2rad(theta0), deg2rad(thetaValues(idx-1)), TR, T1));
    
    if (abs(deg2rad(thetaValues(idx)) - pi/2) < 0.1)
        break;
    end
end

n = idx;

figure(1);
plot(1:n, thetaValues(1:n));
titleText = {'Plot of a variable rf angle as a function of rf pulse number n', ...
    [' for GM at \theta_0 = ', num2str(theta0), '^o '], ...
    [' TR = ', num2str(TR), 'ms and T1 = ', num2str(T1), ' ms']};
title(titleText); 
xlabel('rf pulse number'); ylabel('flip angle (degrees)'); 
grid on




