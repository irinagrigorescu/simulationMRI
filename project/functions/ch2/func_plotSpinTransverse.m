function objHandle = func_plotSpinTransverse (figHandle, vecMu, colors)
% function objHandle = func_plotSpinTransverse (figHandle, vecMu, colors)
%
% a function for plotting a spin
%
% INPUT:
%
% figHandle - the handle to the figure for plotting
%
% vecMu - the magnetic moment to plot
%
% OUTPUT:
%
% objHandle - the handle to the plotted magnetic moment
%
% author: Gary Zhang (gary.zhang@ucl.ac.uk)
%
if nargin < 3
    colors = ['b', 'r'];
end

% Plots the magnetisation vector
objHandle(1) = plot(figHandle, ...
                   [0 vecMu(1)], [0 vecMu(2)], ...
                   'Color', colors(1), 'LineStyle', '-', 'LineWidth', 2);

% Plots the magnetisation vector's tip
objHandle(2) = plot(figHandle, ...
                    vecMu(1), vecMu(2), ...
                    'Color', colors(2), 'Marker', '.', 'MarkerSize', 10);


