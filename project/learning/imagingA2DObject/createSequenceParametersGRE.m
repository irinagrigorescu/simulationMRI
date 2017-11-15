function [dtSampling, Tsampling, dtG, dtPhase, ...
          GxAmplitude, GxArea, ...
          GxPreAmplitude, GxPreArea, ...
          GyAmplitude, dGyAmplitude] = ...
    createSequenceParametersGRE(samplingFrequency, Nx, Ny, ...
                                dkx, dky, gammabar)
% 
% This function creates the kspace parameters
% 
% Input:
%   samplingFrequency = sampling frequency or bandwidth (Hz)
%   Nx,  Ny  = resolution in x and y
%   dkx, dky = dimensions of sampling distances in k-space (1/m)
%   gammabar = reduced gyromagnetic ratio (1/s/T)
% 
% Output:
%   dtSampling = time step between adjacent samples in RO direction (s)
%   Tsampling  = Total sampling duration per line (s)
%   dtG        = duration of half of the readout gradient (s)
%   dtPhase    = duration of phase gradient (s)
% 
%   GxAmplitude    = amplitude of ro gradient (T/m)
%   GxArea         = area of ro gradient (Ts/m)
%   GxPreAmplitude = amplitude of prewinder for ro gradient (T/m)
%   GxPreArea      = area of prewinder for ro gradient (Ts/m)
%   GyAmplitude    = amplitude of first phase encoding gradient (T/m)
%   dGyAmplitude   = gradient step for each phase encoding step (T/m)
% 
% Author: Irina Grigorescu (irina.grigorescu.15@ucl.ac.uk)
% 

% Important time steps
dtSampling = 1/samplingFrequency; 
Tsampling  = dtSampling * Nx; 
dtG        = Tsampling/2;     
dtPhase    = dtG;   

% Readout gradient amplitude and area
GxAmplitude    = dkx/(gammabar*dtSampling);   % T/m
GxArea         = GxAmplitude * Tsampling;     % Ts/m

% Prewinder for RO gradient amplitude and area
GxPreArea      = GxArea / 2;                  % Ts/m
GxPreAmplitude = GxPreArea / dtG;             % T/m

% Phase encoding amplitude and area
dGyAmplitude   = dky/(gammabar*dtPhase);      % T/m
GyAmplitude    = (Ny/2-1)*dGyAmplitude;       % T/m
