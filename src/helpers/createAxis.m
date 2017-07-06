function [xax, yax, zax] = createAxis(xl, yl, zl)
% % % % IRINA GRIGORESCU
% % % % Date created: a long time ago
% % % % Date updated: 03-05-2017
% % % % 
% % % % INPUT:
% % % %     xl, yl, zl - limits of axis of 3D figure
% % % % 
% % % % OUTPUT:
% % % %     xax, yax, zax - the x-/y-/z-axis handles
% % % % 

% % Prerequisites
% Test if all parameters are missing
if nargin < 1
   xl = 1.2;
   yl = 1.2;
   zl = 1.2;
end
% Test if only xl was provided
if nargin < 2
   yl = xl;
   zl = xl;
end
% Test if both xl and yl were provided
if nargin < 3
   zl = min(xl,yl);
end

% Plot z-axis
zax = plot3([[0 0 0]',[0 0 0]'], ...
            [[0 0 0]',[0 0 0]'], ...
            [[5 0 0]',[-5 0 0]'], ...
            '--', 'Color', [0 0 0.1]); % z
hold on
% Plot y-axis
yax = plot3([[0 0 0]',[0 0 0]'], ...
            [[5 0 0]',[-5 0 0]'], ...
            [[0 0 0]',[0 0 0]'], ...
            '--', 'Color', [0 0 0.1]); % y
hold on
% Plot x-axis
xax = plot3([[5 0 0]',[-5 0 0]'], ...
            [[0 0 0]',[0 0 0]'], ...
            [[0 0 0]',[0 0 0]'], ...
            '--', 'Color', [0 0 0.1]); % x
grid on
% Label them
xlabel('x')
ylabel('y')
zlabel('z')
% Limit the axis
xlim([-xl xl])
ylim([-yl yl])
zlim([-zl zl])

end