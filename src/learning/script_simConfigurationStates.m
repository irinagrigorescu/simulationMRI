% % % % IRINA GRIGORESCU
% % % % DATE CREATED: 17-07-2017
% % % % DATE UPDATED: 17-07-2017
% % % % 
% % % % Different configuration states
% % % % 

% % Prerequisites:
clear all
close all
addpath ../helpers/

%% Simulation
% % % % % % % % % % % CREATE SIMULATION:
% % RELAXATION TERMS:
T1 = 1000;
T2 = 100;
% % SPINS:
N   =  5;
% % VIDEO:
videoFlag.flag      = 0;    % % Don't make video
videoFlag.nameVideo = 'configStatesSim-video';

% % % % % % % % % % % RUN SIMULATION:
func_simConfigurationStates(N, T1, T2, videoFlag);

