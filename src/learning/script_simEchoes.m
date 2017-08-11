% % % % IRINA GRIGORESCU
% % % % DATE CREATED: 17-07-2017
% % % % DATE UPDATED: 17-07-2017
% % % % 
% % % % Different Sequences to simulate with the func_simEchoes function
% % % % 

% % Prerequisites:
clear all
close all
addpath ../helpers/

%% 1. CPMG Weigel Paper: 
% % % % % % % % % % % CREATE SIMULATION:
% % RELAXATION TERMS:
T1 =   200;
T2 =    30;
% % EVENTS
rfa = [90, 60, 60, 60, 60, 60, 60, 60]; % flip angles
pha = [90,  0,  0,  0,  0,  0,  0,  0]; % phase angles
tau = [ 0,  2,  4,  4,  4,  4,  4,  4]; % time before them
Ttotal  = 1.2 * sum(tau); % total time of simulation
% % SPINS:
N   =  1000;
dw  =  2.1*pi/(2*pi*min(tau(tau~=0))); % full circle dephasing over 
                                       % smallest tau
dws = linspace(-dw, dw, N); % create a population of spins with equally
                            % distributed off-resonances
% % VIDEO:
videoFlag.flag      = 0;    % % 0 = Don't make video
videoFlag.nameVideo = 'echoSim-CPMG-video';

% % % % % % % % % % % RUN SIMULATION:
func_simEchoes(N, T1, T2, dws, rfa, pha, tau, Ttotal, videoFlag);
                            
%% % % % % % % % 2. Spin Echo Classic
% % RELAXATION TERMS:
T1 =   200;
T2 =    30;
% % % EVENTS
rfa = [90, 150]; % flip angles
pha = [ 0,  0 ]; % phase angles
tau = [ 0,  2 ]; % time before them
Npulses = length(rfa);       % how many
Ttotal  = 4 * sum(tau); % total time of simulation
% % % SPINS:
N   =  1000;
dw  =  2.1*pi/(2*pi*min(tau(tau~=0))); % full circle dephasing over 
                                       % smallest tau
dws = linspace(-dw, dw, N); % create a population of spins with equally
                            % distributed off-resonances
% % VIDEO:
videoFlag.flag      = 0;    % % 0 = Don't make video
videoFlag.nameVideo = 'echoSim-SEclassic-video';

% % % % % % % % % % % RUN SIMULATION:
func_simEchoes(N, T1, T2, dws, rfa, pha, tau, Ttotal, videoFlag);

%% % % % % % % % 3. Spin Echo 8 ball - the 2nd 90 produces a spin echo with
% % reduced amplitude
% % % % EVENTS
% rfa = [90,  90]; % flip angles
% pha = [ 90,   90 ]; % phase angles
% tau = [ 0,   2]; % time before them
% Npulses = length(rfa);       % how many
% Ttotal  = 4 * sum(tau); % total time of simulation
% % % % SPINS:
% N   =  1000;
% dw  =  1.09*pi/(2*pi*min(tau(tau~=0))); % full circle dephasing over 
%                                        % smallest tau
% dws = linspace(-dw, dw, N); % create a population of spins with equally
%                             % distributed off-resonances


%% % % % % % % % 4. Stimulated Echo
% % Magnetization is ?stored? along z-direction 
% % between second and third pulse (so-called ?phase memory?)
% % % % RELAXATION TERMS
% T1 =   200;
% T2 =    30;
% % % % EVENTS
% rfa = [90, 90, 90]; % flip angles
% pha = [ 0,  0,  0 ]; % phase angles
% tau = [ 0,  2, 20]; % time before them
% Npulses = length(rfa);       % how many
% Ttotal  = 2 * sum(tau); % total time of simulation
% % % % SPINS:
% N   =  1000;
% dw  =  1.09*pi/(2*pi*min(tau(tau~=0))); % full circle dephasing over 
%                                        % smallest tau
% dws = linspace(-dw, dw, N); % create a population of spins with equally
%                             % distributed off-resonances


%% % % % % % % % 5. 3 RF PULSES AND 5 ECHOES
% % % % RELAXATION TERMS
% T1 =   200;
% T2 =    80;
% % % % EVENTS
% rfa = [60, 60, 90]; % flip angles
% pha = [ 0,  0,  0 ]; % phase angles
% tau = [ 0,  2,  6]; % time before them
% Npulses = length(rfa);       % how many
% Ttotal  = 2.5 * sum(tau); % total time of simulation
% % % % SPINS:
% N   =  1000;
% dw  =  12.09*pi/(2*pi*min(tau(tau~=0))); % full circle dephasing over 
%                                        % smallest tau
% dws = linspace(-dw, dw, N); % create a population of spins with equally
%                             % distributed off-resonances
                   
% % % % CPMG:
% % rfa = [90, 180, 180, 180, 180, 180]; % flip angles
% % pha = [ 0,   0,   0,   0,   0,   0]; % phase angles
% % tau = [ 0,   5,  10,  10,  10,  10]; % time before them
% % Npulses = length(rfa);       % how many
% % Ttotal  = 1.5 * sum(tau); % total time of simulation

% % % % Carr-Purcell:
% % rfa = [90, 180, 180, 180, 180, 180]; % flip angles
% % pha = [ 0,  90,  90,  90,  90,  90]; % phase angles
% % tau = [ 0,   5,  10,  10,  10,  10]; % time before them
% % Npulses = length(rfa);       % how many
% % Ttotal  = 1.5 * sum(tau); % total time of simulation                            
