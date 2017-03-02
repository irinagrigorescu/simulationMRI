% % % % IRINA GRIGORESCU
% % % % 
% % % % This is a script to test convolution
% % % % 

clear all; close all; clc
addpath ../../helpers

%%
% Properties of the functions
delta_t = .01; % sampling interval
t       = -1:delta_t:1; % time domain spans from -1 to 1
L       = length(t); % number of samples
tp      = 2*t(1):delta_t:2*t(L); % samples for convolution

% Create the functions
ft = zeros(L,1);  gt = zeros(L,1); % the two functions
ft(t>-0.5) = 1;      ft(t>0.5) = 0; % box function
gt(t>-0.5) = t(t>-0.5)+0.5; gt(t>0.5) = 0; % right triangular ramp function

% Calculate convolution and cross-correlation
convfxgx = delta_t .*  conv(ft,gt); % convolution
% Cross-correlation can be written as convolution between the complex
% conjugate of f(-t) and g(t), hence the title in the last plot
corrfxgx = delta_t .* xcorr(ft,gt); % cross-correlation

% % Plot
stepPlot = 4;
colorsPlot = 'rbgk';

figure('Position', [10, 10, 400, 600])
% 1. The functions f(x) and g(x)
subplot(4,1,1)
stem(t(1:stepPlot:end), ft(1:stepPlot:end), colorsPlot(1));
axis([t(1) t(L) 0 2])
xlabel('t'); ylabel('f(t)'); grid on
title({'Looking at convolution and cross-correlation' ; '' ;...
'\quad $f(t) \ast g(t) = \int_{- \infty}^{+ \infty} dx \ f(x) g(t-x)$'}, ...
'Interpreter', 'Latex');

subplot(4,1,2)
stem(t(1:stepPlot:end), gt(1:stepPlot:end), colorsPlot(2));
axis([t(1) t(L) 0 2])
xlabel('t'); ylabel('g(t)'); grid on

% 2. The convolution of f(x) * g(x)
% and the cross-correlation f*(?t) * g(t)
subplot(4,1,3)
stem(tp(1:2*stepPlot:end), convfxgx(1:2*stepPlot:end), colorsPlot(3)),
axis([tp(1) tp(end) 0 1])
xlabel('t'); ylabel({'f(t) \ast g(t)';'(convolution)'});
grid on

subplot(4,1,4)
stem(tp(1:2*stepPlot:end), corrfxgx(1:2*stepPlot:end), colorsPlot(4)),
axis([tp(1) tp(end) 0 1])
xlabel('t'); 
ylabel({'f^*(-t) \ast g(t)';'(cross-correlation)'});
grid on




