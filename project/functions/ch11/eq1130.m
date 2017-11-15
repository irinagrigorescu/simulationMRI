% % % % IRINA GRIGORESCU
% % % % 25/10/2016
% % % % 
% % % % This is a script to understand the Fourier series identity
% % % % Equation 11.30
% % % %

clear all; close all; clc
addpath ../../helpers

%%

% Limits for the series
nLimit = 1000;
% Step for dirac deltas
step = 4;
% Domain of function f(a)
a      = -10:0.1:10;

% Calculate the series
series = 0;
for n = -nLimit:(1/step):nLimit
    series = series + ...
        exp((1i * 2 * pi * n) .* a);
end
% Normalise
series = abs(series) ./ max(abs(series));

% Calculate it as a series of dirac delta functions
% the m will give you the distance between dirac deltas
seriesDirac = 0;
for m = -nLimit:step:nLimit
    seriesDirac = seriesDirac + ...
        dirac_delta_function(a-m, 1e-02);
end

% Plot 
figure('Position', [400, 50, 600, 400])
% Plot the series of Euler formula
stem(a, real(series), 'b*-'); hold on
% Plot the series of dirac delta functions
stem(a, real(seriesDirac), 'gd--');
xlabel('a');
ylabel('f(a)');
title({'$f(a) = \sum_{n=- \infty}^{+ \infty} e^{i 2 \pi n a}$' ; ...
'$\quad \quad = \sum_{m=- \infty}^{+ \infty} \delta (a - m)$'  ; ...
'' }, 'Interpreter', 'Latex');
legend('Euler Formula', 'Dirac Delta', 'Location', 'eastoutside');
xlim([-10 10]);


