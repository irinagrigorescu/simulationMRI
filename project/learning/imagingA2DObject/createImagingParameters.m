function [dkx, dky] = createImagingParameters(Nx, Ny, dx, dy)
% 
% This function creates the kspace parameters
% 
% Input:
%   Nx = resolution of image in x 
%   Ny = resolution of image in y
%   dx = image voxel x dimension (m)
%   dy = image voxel y dimension (m)
% 
% Output:
%   dkx = dimensions of sampling distances in k-space in x direction (1/m)
%   dky = dimensions of sampling distances in k-space in y direction (1/m)
% 
% Author: Irina Grigorescu (irina.grigorescu.15@ucl.ac.uk)
% 

dkx = (2*pi)/(dx*Nx); 
dky = (2*pi)/(dy*Ny);