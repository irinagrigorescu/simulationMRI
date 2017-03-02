% % % % IRINA GRIGORESCU
% % % % DATE: 07-Nov-2016
% % % % 
% % % % CHAPTER 11.2.7
% % % % 
% % % % The derivative theorem
% % % % 

clear all; close all; clc
addpath ../../helpers

%% 
% Plot the 1D Gaussian and its 1st and 2nd derivatives
NSIZE = 1000; % points
xLim  =   20; % domain limit x \in(-xLim, xLim)
x     = linspace(-xLim, xLim, NSIZE); % domain of function
dx    = abs(x(2)-x(3));
dk    = 1/(NSIZE*dx);
k     = linspace(-(NSIZE*dk/2), (NSIZE*dk/2), NSIZE); % k-space domain
%k     = [k linspace(abs(dk), NSIZE*dk/2, NSIZE/2)]; % k-space domain
sigma =  1;    % standard deviation
mu    =  0;    % mean


% The Gaussian function
funcGauss0   = 1/(sigma * sqrt(2*pi)) * ...
    exp(-0.5 .* ((x-mu)./sigma).^2);
% The derivative of the Gaussian function
funcGaussDer = -x./(sigma^3 * sqrt(2*pi)) .* ...
    exp(-0.5 .* ((x-mu)./sigma).^2);
% The second derivative of the Gaussian function
funcGaussDer2 = (sigma^2 - x.^2)./(sigma^5 * sqrt(2*pi)) .* ...
    exp(-0.5 .* ((x-mu)./sigma).^2);


figure
plot(x, funcGauss0)
hold on
plot(x, funcGaussDer)
hold on
plot(x, funcGaussDer2)
grid on
xlabel('$x$', 'Interpreter', 'Latex');
ylabel(['$h(x,\sigma) \ ; \ \frac{\partial h(x,\sigma)}{\partial x} ', ...
    '\ ; \  \frac{\partial h(x,\sigma)}{\partial^2 x}$'], ...
    'Interpreter', 'Latex');
legend('0th order', '1st order', '2nd order');
xlim([-5 5]);

%% 
% Create a function (1D Gaussian) and Fourier transform it
NSIZE = 1000; % points
xLim  =   20; % domain limit x \in(-xLim, xLim)
x     = linspace(-xLim, xLim, NSIZE); % domain of function
dx    = abs(x(2)-x(3));
dk    = 1/(NSIZE*dx);
k     = linspace(-(NSIZE*dk/2), (NSIZE*dk/2), NSIZE); % k-space domain
%k     = [k linspace(abs(dk), NSIZE*dk/2, NSIZE/2)]; % k-space domain
sigma =  1;    % standard deviation
mu    =  0;    % mean


funcGauss  = 1/(sigma * sqrt(2*pi)) * ...
    exp(-0.5 .* ((x-mu)./sigma).^2);
funcGaussT = exp((-(k.^2).*(sigma^2)) ./ 2);

% e^f(x)' = f'(x) e^f(x)
funcGaussDeriv = -x./(sigma^3 * sqrt(2*pi)) .* ...
    exp(-0.5 .* ((x-mu)./sigma).^2);

funcGaussDeriv2 = (sigma^2 - x.^2)./(sigma^5 * sqrt(2*pi)) .* ...
    exp(-0.5 .* ((x-mu)./sigma).^2);


figure
subplot(3,2,1)
plot(x, funcGauss)
grid on
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$h(x,\sigma)$', 'Interpreter', 'Latex')
title('$\frac{1}{\sigma \sqrt{2 \pi}} e^{- \frac{x^2}{2 \sigma^2} }$', 'Interpreter', 'Latex');
xlim([-5 5])


subplot(3,2,2)
plot(k, funcGaussT)
xlim([-5, 5])
grid on
xlabel('$k$', 'Interpreter', 'Latex')
ylabel('$F[h(x,\sigma)]$', 'Interpreter', 'Latex')
title('$H(k) = e^{-\frac{k^2 \sigma^2}{2}}$', 'Interpreter', 'Latex');


subplot(3,2,3)
plot(x, funcGaussDeriv)
grid on
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$\frac{\partial h(x,\sigma)}{\partial x}$', 'Interpreter', 'Latex');
title('$-\frac{x}{\sigma^3 \sqrt{2 \pi}} e^{- \frac{x^2}{2 \sigma^2} }$', 'Interpreter', 'Latex')
xlim([-5 5])


subplot(3,2,4)
plot(k, abs(1i*2*pi.*k.*funcGaussT))
hold on
plot(k, real(1i*2*pi.*k.*funcGaussT), '--')
hold on
plot(k, imag(1i*2*pi.*k.*funcGaussT), '--')
xlim([-5 5])
grid on
xlabel('$k$', 'Interpreter', 'Latex')
ylabel('$F[\frac{\partial h(x,\sigma)}{\partial x}]$', 'Interpreter', 'Latex')
title('$i 2 \pi k H(k)$', 'Interpreter', 'Latex')
legend('abs','real','imag')

subplot(3,2,5)
plot(x, funcGaussDeriv2)
grid on
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$\frac{\partial h(x,\sigma)}{\partial^2 x}$', 'Interpreter', 'Latex');
title('$-\frac{\sigma^2 - x^2}{\sigma^5 \sqrt{2 \pi}} e^{- \frac{x^2}{2 \sigma^2} }$', 'Interpreter', 'Latex')
xlim([-5 5])


subplot(3,2,6)
plot(k, abs(((1i*2*pi.*k).^2).*funcGaussT))
hold on
plot(k, real(((1i*2*pi.*k).^2).*funcGaussT), '--')
hold on
plot(k, imag(((1i*2*pi.*k).^2).*funcGaussT), '--')
xlim([-5 5])
grid on
xlabel('$k$', 'Interpreter', 'Latex')
ylabel('$F[\frac{\partial h(x,\sigma)}{\partial^2 x}]$', 'Interpreter', 'Latex')
title('$(i 2 \pi k)^2 H(k)$', 'Interpreter', 'Latex')
legend('abs','real','imag')
