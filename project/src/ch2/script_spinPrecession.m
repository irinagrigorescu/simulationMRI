% 
% This script shows spin precession
% 
% Author: Irina Grigorescu, irina.grigorescu.15@ucl.ac.uk
%                           irinagry@gmail.com
% 
% Function dependencies:
%   plotVectorComponents
%   plotVectorFromOrigin
%   Rotz
%   func_saveToVideo

clear all ; close all ; clc
addpath(genpath('../../helpers/'))

%% Prerequisites regarding saving your animation as a video mp4 file
videoFlag = struct;
videoFlag.flag = 0;
videoFlag.nameVideo = '../../figures/videos/spinPrecessionAnimation';

% Create a video writer class object 
videoClass = VideoWriter([videoFlag.nameVideo, '.mp4'], ...
                         'MPEG-4');    % #need this for video
frameRate = 2;                         % optional
videoClass.set('FrameRate',frameRate); % optional
open(videoClass)                       % #need this for video

%% The initial position of the magnetic moment vector 
% Polar coordinates of the position:
theta = pi/4; % polar angle
psi   = pi/2; % azimuthal angle
% Magnitude of vector:
muMagnitude = 1;
% Components of vector:
muVec = muMagnitude * [ cos(psi) * sin(theta) ; ...
                        sin(psi) * cos(theta) ; ...
                              cos(theta)      ] ;

%% The precession information
% Precession frequency (rad*kHz) (relative to Larmor)
dOmega = 2*pi;
% Time increment (ms)
dT = 0.01;
% Angle of rotation for each time increment
dAngle = dOmega * dT;
% Number of increments
nIncrements = 500;

%% Calculate the magnetisation vector position for each time step
muVecHistory = zeros(3, nIncrements+1);
muVecHistory(:,1) = muVec;

for i = 2:nIncrements+1
    muVecHistory(:, i) = Rotz(-dAngle) * muVecHistory(:, i-1); 
end

%% Plot function as an animation for each time step
figHandle = figure('Position', [50, 150, 1200, 600]); hold on

toClose = 0;

for i = 1:nIncrements+1
    % Plot the 3D vector animation
    if i == 1
        fh1 = subplot(3,4,[1, 10]);
        % Create the 3D axis
        createAxis(muMagnitude); view(135, 15); axis square
        % title
        title('Animation of spin precession')
        % axis labels
        xlabel('\mu_x'); ylabel('\mu_y'); zlabel('\mu_z');
    end
	% Plot vector
    [handleLine, handleTip] = plotVectorFromOrigin(fh1, muVecHistory(:,i));
	
    
    % Plot the component-wise animation
    if i == 1
        fh2 = subplot(3,4,[7,8]);
        % plot component wise
        plotVectorComponents(fh2, i*dT, muVecHistory(:,i))
        % legend
        legend('\mu_x','\mu_y','\mu_z', 'Location', 'eastoutside')
        % plot limits
        xlim([0 nIncrements*dT]); 
        ylim([-muMagnitude muMagnitude]);
        % title
        title('The magnetic moment vector component-wise')
        % plot labels
        xlabel('time (ms)'); ylabel('\mu_i');
    end
    % plot component wise
    plotVectorComponents(fh2, i*dT, muVecHistory(:,i))
    
    % Send current figure handle to the saveToVideo function to be saved
	func_saveToVideo(videoClass, figHandle, videoFlag, toClose); % set to 0 
                                                               % to not yet 
                                                               % close the 
                                                             % video handle
    
    % Pause
    pause(0.01)
    
    % Change color for historical reasons
    if i ~= (nIncrements+1)
        handleTip.Color(1) = 0.9;
        handleTip.Color(2) = 0.9;
        handleTip.Color(3) = 0.9;
        handleTip.MarkerSize = 10;
        handleLine.Color(4) = 0.1;
    end
    
    % Next iteration close the video handle
    if i == nIncrements
        toClose = 1;
    end
   
end
