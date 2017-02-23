% % % % IRINA GRIGORESCU
% % % % DATE: 13-Feb-2017
% % % % CHAPTER 18.1.1
% % % % 
% % % % This script plots the spin density behaviour 
% % % % at TE = 0 for rf flip angles ranging from 0 to pi
% % % % with TR = 40ms and for WM,GM,CSF,Fat at 1.5T
% % % % 

%% Loading some data
addpath ../../helpers
run config; params = ans; clear ans;

%% Signal function 
% This is the steady-state signal from 
% a voxel containing one isochromat (T2 relaxation)
% or several isochromats (T2* relaxation)
% This is eq 18.4 Green Bible
rhohat = @(theta, TE, TR, rho0, T1, T2) ...
    rho0 .* sin(theta) .* ...
    (( ( 1 - exp(-TR/T1) ) ./ (1 - cos(theta).*exp(-TR/T1)) ) * ...
    exp(- TE / T2));

%% Plot figure
% Plotting the signal for 0-pi range of rf pulse angles
% for 4 tissue types (csf,wm,gm,fat)

% RF pulse angles ranging from 0 to pi
theta = 0:0.01:pi/4; 
% A range of TEs
te = 0:0.1:40;
% Learning how many tissue types my structure has
tissueNo = length(fieldnames(params.tissue));
% The figure legend is constructed here
legendMat = cell(tissueNo,1);
% Keep signal values in a structure of tissues 
% e.g. if 4 tissue types then 
% rhoTissue = structure of 4 same-length vectors of signal values
rhoTissueByTheta = struct;
% Similar to above but for set angle and a range of TEs
rhoTissueByTEs = struct;


% For each tissue type
for tissueIdx = 1:tissueNo
    
    % Get Names of tissues to refer to them later
    tissueType = fieldnames(params.tissue);
    
    % Get values for plot for each tissue type
    % for TE = 0ms and TR = 40ms
    rhoTissueByTheta.(tissueType{tissueIdx}) = ...
        rhohat(theta, 0, 40, ...
        params.tissue.(tissueType{tissueIdx}).rho, ...
        params.tissue.(tissueType{tissueIdx}).t1, ...
        params.tissue.(tissueType{tissueIdx}).t2);
    
    figure(1)
    plot(theta * 180/pi, rhoTissueByTheta.(tissueType{tissueIdx}), '-');
    hold on

    % Get values for plot for each tissue type
    % for a range of TEs, theta = pi/2 and TR = 5000ms
    rhoTissueByTEs.(tissueType{tissueIdx}) = ...
        rhohat(pi/10.58, te, 40, ...
        params.tissue.(tissueType{tissueIdx}).rho, ...
        params.tissue.(tissueType{tissueIdx}).t1, ...
        params.tissue.(tissueType{tissueIdx}).t2);
    
    figure(2)
    plot(te, rhoTissueByTEs.(tissueType{tissueIdx}), '-');
    hold on
    
    legendMat{tissueIdx} = tissueType{tissueIdx};
end

figure(1); 
legend(legendMat); 
title('Comparison of the SSI signal for different RF angles and TE = 0ms'); 
xlabel('\theta'); ylabel('\rho(\theta,TE)'); 

figure(2); 
legend(legendMat); 
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
        plot(angleOfIntersect, intersect, 'b*')
    end
    
end








