% % % % IRINA GRIGORESCU
% % % % DATE: 08-Nov-2016
% % % %
% % % % CHAPTER 9.28
% % % % 
% % % % The Dirac Delta function can be written as 
% % % % lim K->infty [I(x,K)] = lim K->infty [2K sinc(2 pi K x)]

clear all; close all; clc;
addpath ../../helpers

%% 
% Create the I(z,K) function
NSIZE = 10000;   % points
xLim  =     0.2; % domain limits x \in(-xLim, xLim)
x     = linspace(-xLim, xLim, NSIZE); % domain of function
K = [1 5 10 50 100]; % a few values for the K constant


% Create the I(x,K) function 
funcDiracI = zeros(numel(K), NSIZE); % for each constant K

for i = 1:numel(K)
    funcDiracI(i,:) = 2*K(i) .* sinc((2 * pi * K(i)) .* x);
end

%%
% Plot the function for different values of K
% It can be seen that for K approaching infinity the function
% approximates the Dirac Delta increasingly well

legendText = cell(1,numel(K)); % legend text

% Plot the function for different K values
for i = 1:numel(K)
    figure(1)
    plot(x, (funcDiracI(i,:)));
    legendText{i} = sprintf('K = %d', K(i));
    hold on
end
legend(legendText);
xlabel('x',      'Interpreter','Latex');
ylabel('I(x,K)', 'Interpreter','Latex');
title('$I(x,K) = 2K sinc(2 \pi K x)$', 'Interpreter', 'Latex')
xlim([-xLim/5 xLim/5]);



