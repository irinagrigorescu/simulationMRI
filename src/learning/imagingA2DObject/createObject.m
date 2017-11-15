function [object2D, mask] = createObject(FOVx, FOVy, T1, T2, flagPlot)
% 
% This function creates an object with position in x, y and material
% properties described by T1 and T2
% 
% Input:
%   FOVx     = size of object in x (mm)
%   FOVy     = size of object in y (mm)
%   T1       = array of T1s (s)
%   T2       = array of T2s (s)
%   flagPlot = 1|0 to plot or not to plot
% 
% Output:
%   object2D = object represented by FOVx x FOVy x 4
%              4 = xcoord (mm), ycoord (mm), T1 (s), T2 (s)
%   mask = structure with FOVx x FOVy masks 
% 
% Author: Irina Grigorescu (irina.grigorescu.15@ucl.ac.uk)
% 

% Create positions 
posx = (-FOVx/2)/1E03 : 1/1E03 : (FOVx/2-1)/1E03; % x positions in meters
posy = (-FOVy/2)/1E03 : 1/1E03 : (FOVy/2-1)/1E03; % y positions in meters
[Xcoord, Ycoord] = meshgrid(posx, posy); Ycoord = Ycoord.'; % coordinates

% My object has 4 parameters: xcoord, ycoord, T1, T2
object2D = zeros(FOVx, FOVy, 4); 
object2D(:,:,1) = Xcoord;
object2D(:,:,2) = Ycoord;

% Create masks for 3 objects in the image
mask1 = createCirclesMask([FOVx FOVy], ...
                          [FOVx/2-FOVx/4, FOVy/2-FOVy/4], 8);
                                                 % upper left circle
mask2 = createCirclesMask([FOVx FOVy], ...
                          [FOVx/2, FOVy/2], 8); % centre circle
mask3 = createCirclesMask([FOVx FOVy], ...
                          [FOVx/2+FOVx/4, FOVy/2+FOVy/4], 8);
                                                 % lower right circle

% Create a structure that keeps all the masks together
mask = struct;
mask.mask1 = mask1;
mask.mask2 = mask2;
mask.mask3 = mask3;
                                                 
% Create object
object2D(:,:,3) = mask1.*T1(4);
object2D(:,:,4) = mask1.*T2(4);

% Figure:
if flagPlot == 1
    figure
    subplot(1,2,1)
    imagesc(squeeze(object2D(:,:,3)));
    axis off; axis square; colorbar
    title('T_1')

    subplot(1,2,2)
    imagesc(squeeze(object2D(:,:,4)))
    axis off; axis square; colorbar
    title('T_2')
end
