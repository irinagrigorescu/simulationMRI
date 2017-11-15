% % % % IRINA GRIGORESCU
% % % % DATE: 10-Nov-2016
% % % % 
% % % % CHAPTER 11.2.8
% % % % 
% % % % Fourier Transform Symmetries
% % % % 

clear all; close all; clc
addpath ../../helpers

%% 
% Plot the sine and cosine functions
NSIZE = 1001; % points
xLim  =   20; % domain limit x \in(-xLim, xLim)
x     = linspace(-xLim, xLim, NSIZE); % domain of function
dx    = abs(x(2)-x(3));
dk    = 1/(NSIZE*dx);
k     = linspace(-(NSIZE*dk/2), NSIZE*dk/2, NSIZE); % k-space domain
%k     = [k linspace(abs(dk), NSIZE*dk/2, NSIZE/2)]; % k-space domain

k_0 = 1;

sinFunc = sin(2*pi*k_0.*x);
cosFunc = cos(2*pi*k_0.*x);


%% Plotting the sine and cosine functions and their fourier transforms
figure, 

subplot(2,2,1)
plot(x, cosFunc)
title('$cos(2 \pi k_0 x)$', 'Interpreter', 'Latex')
xlabel('x')
grid on
xlim([-2 2])

subplot(2,2,2)
plot(k, myNorm(real(fftshift(fft(cosFunc)))))
hold on
plot(k, myNorm(imag(fftshift(fft(cosFunc)))), '--')
title(['\textit{F}$[cos(2 \pi k_0 x)] = \frac{1}{2} \delta(k-k_0) ', ...
    '+ \frac{1}{2} \delta(k+k_0)$'], 'Interpreter', 'Latex');
xlabel('k')
grid on
xlim([-2 2])

subplot(2,2,3)
plot(x, sinFunc)
title('$sin(2 \pi k_0 x)$', 'Interpreter', 'Latex')
xlabel('x')
grid on
xlim([-2 2])

subplot(2,2,4)
plot(k, myNorm(real(fftshift(fft(sinFunc)))))
hold on
plot(k, myNorm(imag(fftshift(fft(sinFunc)))), '--')
title(['$Im\{F[sin(2 \pi k_0 x)]\} = \frac{1}{2} \delta(k+k_0)',...
    ' - \frac{1}{2} \delta(k-k_0)$'], 'Interpreter', 'Latex');
xlabel('k')
grid on
xlim([-2 2])






