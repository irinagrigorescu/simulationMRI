function plotVectorComponents(figureHandle, tp, vectorComponents)
% This function plots the vector components 
%
% INPUT:
%   figureHandle = the figure handle in which to plot
%   tp = time point at which to plot
%   vectorComponents = the vector you wish to plot
% OUTPUT:
%   -
%
% Author: Irina Grigorescurina.grigorescu.15@ucl.ac.uk
%                           irinagry@gmail.com

colors = ['r', 'g', 'b', 'k'];

for i = 1:3
	plot(figureHandle, ...
         tp, vectorComponents(i), ...
         'Color', colors(i), 'Marker', '.', 'MarkerSize', 2);
    hold on
end
