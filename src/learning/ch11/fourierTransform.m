% % % % Playing with the Fourier transform
% % % % 
% % % % f(t) <=>      F(w)
% % % % 

clear all; close all; clc
addpath ../../helpers

%%

% Parameters for my function in the time domain
N  = 1001; % number of points

Td =   5; % sampling time (from 0 to Td)
x  = linspace(0, Td, N-1); % domain of function
fx = cos(2*pi*x) + 4.*cos(7*2*pi*x);
dx = Td/N;

Fs  = N/Td; %1/dx;    % sampling frequency
df  = 1/Td;           % df
xF  = 0:df:Fs/2;  % domain of function
% xF  = linspace(0, Fs/2, N/2+1);
fxF = fftshift(fft(fx)); % function fourier transformed



%%
figure, 

% Plot f(x)
subplot(2,1,1)
plot(x, fx)
xlabel('x')
ylabel('f(x)')
grid on

% Plot F[f(x)]
subplot(2,1,2)
FtoShow = 10/df; % Show up until frequnecy 10
plot(xF(1:FtoShow+1), ...
     abs(fxF(round(N/2) : round(N/2) + FtoShow) ./ N)) 
        % normalised by number of samples
xlabel('1/x')
ylabel('F[f(x)]')
grid on


%% Fourier transforms on some functions
% Amplitudes 
A =  8; 
B = 10;
% Number of sampling points 
L = 500;
% Domain of function
t = linspace(0,1,L);
mask = ones(1,L); mask(2:50:end) = 0;
% Functions
fcos  = A * cos(2*pi*5*t); %+ B * cos(2*pi*10*t) + A2 * cos(2*pi*5*t);
fsin = A * sin(2*pi*5*t);
fm = A * cos(2*pi*6*t) + B * cos(2*pi*10*t);
fmix = 0; N = 10;
for i = 1:N
    fmix = fmix + cos(2*pi*(2*i)*t);
end

% Fourier transform the function
fcosF  = fftshift(fft(fcos));
fsinF  = fftshift(fft(fsin));
fmF    = fftshift(fft(fm));
fmixF  = fftshift(fft(fmix));

% % % % PLOT
% Plot time domain and frequency domain
figure('Position', [50,50,1000,600] ),
plotL = 3;
plotC = 2;

% % % % % % % % % % % % 
% Time domain functions

% % Sine and Cosine
subplot(plotL,plotC,1)
plot(t,fcos)
hold on
plot(t,fsin)
axis([0, 1, -9, 9]), grid on
%title([sprintf('Time Domain: \n'), ...
title([...
      '$$cos(2 \pi s t)$$', ' and ', '$$sin(2 \pi s t)$$'], ...
      'interpreter', 'latex');
legend([num2str(A), ' cos(2 \pi 5 t)' ],...
       [num2str(A), ' sin(2 \pi 5 t)'], ...
       'Location', 'northeast');

% % Mixture of sine and cosine
subplot(plotL,plotC,3)
plot(t,fm)
axis([0, 1, -20, 20]), grid on
title('$$A cos(2 \pi s_1 t) + B cos(2 \pi s_2 t)$$', ...
      'interpreter', 'latex');
legend([num2str(A), ' cos(2 \pi 6 t) + ', ...
        num2str(B), ' cos(2 \pi 10 t)'], ...
       'Location', 'northeast');
   
% % Mixture of increasing frequencies
subplot(plotL,plotC,5)
plot(t,fmix)
axis([0, 1, -5, 11]), grid on
title(['$$\sum_{k=1}^{k=',num2str(N),'} cos(2\pi (2k) t)$$'], ...
      'interpreter', 'latex');
   
   
% % % % % % % % % % % % % % 
% Frequency domain function

% % Sine and Cosine
subplot(plotL,plotC,2)
plot((-L/2 : L/2-1), abs(fcosF/L), '-')
hold on
plot((-L/2 : L/2-1), abs(fsinF/L), '--')
axis([-L/20, L/20, -2, 5]), grid on
title([...
      '$$F(cos(2 \pi s t))$$', ' and ', '$$F(sin(2 \pi s t))$$'], ...
      'interpreter', 'latex');
   
% % Mixture of sine and cosine
subplot(plotL,plotC,4)
plot((-L/2 : L/2-1), abs(fmF/L), '-')
axis([-L/20, L/20, -2, 7]), grid on
title('$$F(A cos(2 \pi s_1 t) + B cos(2 \pi s_2 t))$$', ...
      'interpreter', 'latex');

% % Mixture of increasing frequencies
subplot(plotL,plotC,6)
plot((-L/2 : L/2-1), abs(fmixF/L), '-')
axis([-L/20, L/20, -.5, 1]), grid on
title(['$$F( \sum_{k=1}^{k=',num2str(N),'} cos(2\pi (2k) t))$$'], ...
      'interpreter', 'latex');


