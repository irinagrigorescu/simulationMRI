% % % % IRINA GRIGORESCU
% % % % DATE: 15-Feb-2017
% % % % CHAPTER 18.1.3 Approach to incohere steady state
% % % % 
% % % % This script plots the normalised longitudinal magnetization
% % % % of four tissue types as a function of the rf pulse number
% % % % for 2 different TRs and a fixed flip angle
% % % % 

clear all; close all; clc
addpath ../../helpers

%% Loading some data
run config; params = ans; clear ans;

%% The relaxation function
E12 = @(tn,t12) exp(-tn./t12);

%% Function that calculates the steady-state or equilibrium value
%  Mze for a given rf flip angle and tissue type (T1)
Mze = @(tn,rho0,t1,thetaAngle) (rho0 .* (1 - E12(tn, t1)))          ./ ...
                               (1 - E12(tn, t1) .* cos(thetaAngle));

%% Function that calculates the longitudinal magnetization behaviour
%  as a function of rf flip angle and the rf pulse cycle
Mz = @(n,tr,rho0,t1,thetaAngle) ...
                (rho0 - Mze(tr,rho0,t1,thetaAngle)) .* ((E12(tr,t1)).^n) + ...
                 Mze(tr,rho0,t1,thetaAngle) ;  
                   
%% Plot figure
TR1  =  40;
TR2  = 400;
N    = 200;
theta = 10;
% Learning how many tissue types my structure has
tissueNo = length(fieldnames(params.tissue));
% The figure legend is constructed here
legendMat = cell(tissueNo,1);

% For each tissue type
for tissueIdx = 1:tissueNo
    % Get Names of tissues to refer to them later
    tissueType = fieldnames(params.tissue);
    
    % Get values for plot for each tissue type
    % for TR = 40ms
    MzByTheta.(tissueType{tissueIdx}) = ...
        Mz((1:N), TR1, ...
        params.tissue.(tissueType{tissueIdx}).rho, ...
        params.tissue.(tissueType{tissueIdx}).t1, ...
        theta*pi/180);
    
    figure(1)
    subplot(1,2,1)
    plot((1:N), MzByTheta.(tissueType{tissueIdx}), '-');
    hold on
    
    % Get values for plot for each tissue type
    % for TR = 400ms
    MzByTheta.(tissueType{tissueIdx}) = ...
        Mz((1:N), TR2, ...
        params.tissue.(tissueType{tissueIdx}).rho, ...
        params.tissue.(tissueType{tissueIdx}).t1, ...
        theta*pi/180);
    
    subplot(1,2,2)
    plot((1:N), MzByTheta.(tissueType{tissueIdx}), '-');
    hold on
   
    legendMat{tissueIdx} = tissueType{tissueIdx};
end

figure(1);
subplot(1,2,1)
legend(legendMat); 
titleText = {'Plot of M_z^-(n) as a function of rf pulse number n', ...
    [' for different tissue types at \theta = ', num2str(theta), '^o '], ...
    [' and TR = ', num2str(TR1), 'ms']};
title(titleText); 
xlabel('n'); ylabel('normalised M_z^-(n)'); 
grid on


figure(1);
subplot(1,2,2)
legend(legendMat); 
titleText = {'Plot of M_z(n) as a function of rf pulse number n', ...
    [' for different tissue types at \theta = ', num2str(theta), '^o '], ...
    [' and TR = ', num2str(TR2), 'ms']};
title(titleText); 
xlabel('n'); ylabel('normalised M_z(n)');
grid on

