% % % % IRINA GRIGORESCU
% % % % DATE: 13-Feb-2017
% % % % CHAPTER 18.1.1
% % % % 
% % % % This script plots the spin density behaviour 
% % % % at TE = 0 for rf flip angles ranging from 0 to pi
% % % % with TR = 40ms and for WM,GM,CSF,Fat at 1.5T
% % % % 

clear all; close all; clc;
addpath ../../helpers

%% Loading some data
run config; params = ans; clear ans;

%% Functions
% % Signal function 
% This is the steady-state signal from 
% a voxel containing one isochromat (T2 relaxation)
% or several isochromats (T2* relaxation)
% This is eq 18.14 Green Bible
rhohatFunc = @(theta, TE, TR, rho0, T1, T2) ...
    rho0 .* sin(theta) .* ...
    (( ( 1 - exp(-TR/T1) ) ./ (1 - cos(theta).*exp(-TR/T1)) ) * ...
    exp(- TE / T2));

% % Longitudinal magnetisation evolution
mzeFunc = @(theta, TE, TR, rho0, T1, T2) ...
    1 .* ...
    ( ( 1 - exp(-TR/T1) ) ./ (1 - cos(theta).*exp(-TR/T1)) );

%% Plot figure
% Plotting the signal for 0-pi range of flip angles
% for 4 tissue types (csf,wm,gm,fat)

% Flip angles ranging from 0 to pi
theta = 0:0.01:pi/4; 
% A range of TEs
te = 0:0.1:40;
% Learning how many tissue types my structure has
tissueNo = length(fieldnames(params.tissue));
% The figure legend is constructed here
legendMat = cell(tissueNo*2,1);
% Keep signal values in a structure of tissues 
% e.g. if 4 tissue types then 
% rhoTissue = structure of 4 same-length vectors of signal values
rhoTissueByTheta = struct;
% Keep longitudinal magnetization values in a structure of tissues 
% e.g. if 4 tissue types then 
% rhoTissue = structure of 4 same-length vectors of signal values
mzeTissueByTheta = struct;
% Similar to above but for set angle and a range of TEs
rhoTissueByTEs = struct;

plotColors = 'rgbk';

% For each tissue type
for tissueIdx = 1:tissueNo
    
    % Get Names of tissues to refer to them later
    tissueType = fieldnames(params.tissue);
    
    % Get values for plot for each tissue type
    % for TE = 0ms and TR = 40ms
    rhoTissueByTheta.(tissueType{tissueIdx}) = ...
        rhohatFunc(theta, 0, 40, ...
            params.tissue.(tissueType{tissueIdx}).rho, ...
            params.tissue.(tissueType{tissueIdx}).t1, ...
            params.tissue.(tissueType{tissueIdx}).t2);
    
    mzeTissueByTheta.(tissueType{tissueIdx}) = ...
        mzeFunc(theta, 0, 40, ...
            params.tissue.(tissueType{tissueIdx}).rho, ...
            params.tissue.(tissueType{tissueIdx}).t1, ...
            params.tissue.(tissueType{tissueIdx}).t2);
    
    % Plot signal by tissue type
    figure(1)
    %subplot(1,2,1)
    plot(rad2deg(theta), rhoTissueByTheta.(tissueType{tissueIdx}), ...
                ['-', plotColors(tissueIdx)]);
    hold on

    % Plot longitudinal magnetisation by tissue type
    figure(1)
    %subplot(1,2,2)
    plot(rad2deg(theta), mzeTissueByTheta.(tissueType{tissueIdx}), ...
                ['--', plotColors(tissueIdx)]);
    hold on
    
    % Get values for plot for each tissue type
    % for a range of TEs, theta = pi/2 and TR = 5000ms
    rhoTissueByTEs.(tissueType{tissueIdx}) = ...
        rhohatFunc(pi/10.58, te, 40, ...
        params.tissue.(tissueType{tissueIdx}).rho, ...
        params.tissue.(tissueType{tissueIdx}).t1, ...
        params.tissue.(tissueType{tissueIdx}).t2);
    
    % Plot signal values for Ernst angle
    figure(2)
    plot(te, rhoTissueByTEs.(tissueType{tissueIdx}), ...
                ['-', plotColors(tissueIdx)]);
    hold on
    
    legendMat{(tissueIdx-1)*2+1} = tissueType{tissueIdx};
    legendMat{(tissueIdx-1)*2+2} = tissueType{tissueIdx};
end

figure(1); 
%subplot(1,2,1)
legend(legendMat); 
title('Comparison of the SSI signal for different flip angles and TE = 0ms'); 
xlabel('\theta'); ylabel('\rho(\theta,TE) and M_z_e(\theta,TE)'); 
% %subplot(1,2,2)
% legend(legendMat); 
% title('Comparison of the SSI signal for different flip angles and TE = 0ms'); 
% xlabel('\theta'); ylabel('\rho(\theta,TE)'); 

figure(2); 
legend(legendMat{1:2:end}); 
title('Comparison of the SSI signal for different TEs and fixed \theta = 17^o');
xlabel('TE (ms)'); ylabel('\rho(\theta,TE)');

%% Find intersection of WM and GM curves
for i = 2:length(theta)
    pair1 = rhoTissueByTheta.(tissueType{2})(i-1) - rhoTissueByTheta.(tissueType{3})(i-1);
    pair2 = rhoTissueByTheta.(tissueType{2})(i)   - rhoTissueByTheta.(tissueType{3})(i);
       
    if sum(sign([pair1 pair2])) == 0
        angleOfIntersect = mean([theta(i), theta(i-1)]) .* 180/pi;
        intersect = abs(rhoTissueByTheta.(tissueType{2})(i-1) + rhoTissueByTheta.(tissueType{2})(i) ) / 2;
        figure(1)
        hold on
        %subplot(1,2,1)
        plot(angleOfIntersect, intersect, 'b*')
    end
    
end








