% % % % IRINA GRIGORESCU
% % % % DATE: 07-Nov-2016
% % % % 
% % % % CHAPTER 11.2.4
% % % % 
% % % % The Duality theorem 
% % % % 
% % % % Duality between the time and frequency domains is another important
% % % % property of Fourier transforms. 
% % % % 
% % % % f(t) <=>      F(w)
% % % % F(t) <=> 2 pi f(-w)
% % % % 

clear all; close all; clc
addpath ../../helpers

%%
% % % % This script is to see the periodicity of H(k)

% Parameters for my function
NSIZE  = 1000; % number of points
kLimit =   2; % domain limit k \in(-kLimit, kLimit)
k   = linspace(-kLimit, kLimit, NSIZE); % domain of function
x_0 = [0 1 2]; % a few values for the x_0 constant

%
% Create function H(k) = e^(-i 2pi k x_0)
Hk = zeros(numel(x_0), NSIZE); % for each constant x_0

for i = 1:numel(x_0)
    Hk(i,:) = exp(- 1j .* 2 .* pi .* k .* x_0(i));
end

%
% Plot the function for different x_0 values
textLegend = cell(1,numel(x_0)); % legend text
for i = 1:numel(x_0)
    figure(1)
    plot(k, real(Hk(i,:)));
    hold on;
    textLegend{i} = ['x_0 = ', num2str(x_0(i))];
end
grid on

% % Adding text to the figure
text(-0.50, 0.9,'|\leftarrow', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 14);
text(-0.25, 0.9,'1/x_0', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 12);
text(-0.25, 0.8,'(for x_0 = 2)', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 10);
text( 0.00, 0.9,'\rightarrow|', ...
    'HorizontalAlignment', 'right', ...
    'FontSize', 14);
xlabel('k');
ylabel('H(k)');
title('$H(k) = e^{- i 2 \pi k x_0}$', 'Interpreter', 'Latex');
legend(textLegend);


%%
% % % % Create the famous sinc(x) / rect(?) fourier transform pair.
% W sinc(pi W k) <-> rect(x/W)


% Parameters for my functions
% number of points
NSIZE = 10000; 
% domain limit x \in(-xLimit, xLimit)
xLim  =     5; 
% domain of function
x     = linspace(-xLim, xLim, NSIZE); 
dx    = abs(x(2)-x(3));
dk    = 1/(NSIZE*dx);
% domain of fourier transf
k     = linspace(-(NSIZE*dk/2), (NSIZE*dk/2), NSIZE);

% Create the functions
W         = 2; % for sinc function
sincFunc  = W.*sinc(pi*W*k);
sincFuncT = W.*sinc(pi*W.*x);
rectFunc  = zeros(NSIZE,1);
rectFuncT = zeros(NSIZE,1);
for i = 1:NSIZE
    if x(i) > -W && x(i) < W
        rectFunc(i) = 1;
    end
    if k(i) > -W && k(i) < W
        rectFuncT(i) = 1;
    end
end

figure
subplot(2,2,1)
plot(x, rectFunc)
grid on
xlabel('x')
ylabel('$rect(\frac{x}{W})$', 'Interpreter', 'Latex')


subplot(2,2,2)
plot(k, sincFunc); %real(fftshift(fft(rectFunc))./NSIZE))
xlim([-5, 5])
grid on
xlabel('k')
ylabel('$W \ sinc(\pi W k)$', 'Interpreter', 'Latex')


subplot(2,2,3)
plot(x, sincFuncT)
grid on
xlabel('x')
ylabel('$W \ sinc(\pi W x)$', 'Interpreter', 'Latex')

subplot(2,2,4)
plot(k, rectFuncT) %abs(fftshift(fft(W.*sinc(pi*W.*x)))./NSIZE))
xlim([-5, 5])
grid on
xlabel('k')
ylabel('$rect(\frac{k}{W})$', 'Interpreter', 'Latex')

















