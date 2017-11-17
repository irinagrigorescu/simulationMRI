%% Function that plots the Figure3.2

addpath ../helpers

%createAxis();
view([-45 -45 10])

% Draw unprimed axis
quiver3(0,0,0,0,0,1, 'Color', [0,0,0]), hold on
quiver3(0,0,0,0,-1,0, 'Color', [0,0,0]), hold on
quiver3(0,0,0,-1,0,0, 'Color', [0,0,0]), hold on

% Draw text
text(0.1, 0.1, 0.7, '$z$', 'FontSize', 24, 'Interpreter', 'latex'), hold on
text(0.1, -0.7, 0.1, '$y$', 'FontSize', 24, 'Interpreter', 'latex'), hold on
text(-0.7, 0.1, 0.1, '$x$', 'FontSize', 24, 'Interpreter', 'latex'), hold on

% Draw primed axis
quiver3(0,0,0,-1,0,1, 'Color', [0.7,0,0.7]), hold on
quiver3(0,0,0,0,-1,1, 'Color', [0.7,0,0.7]), hold on
quiver3(0,0,0,-1,-1,0, 'Color', [0.7,0,0.7]), hold on

% Draw text
text(-0.7, 0.1, 0.7, '$z''$', 'FontSize', 24, 'Color', [0.7,0,0.7], ...
    'Interpreter', 'latex'), hold on
text(0.1, -0.9, 0.8, '$y''$', 'FontSize', 24, 'Color', [0.7,0,0.7], ...
    'Interpreter', 'latex'), hold on
text(-0.8, -0.9, 0.1, '$x''$', 'FontSize', 24, 'Color', [0.7,0,0.7], ...
    'Interpreter', 'latex'), hold on


