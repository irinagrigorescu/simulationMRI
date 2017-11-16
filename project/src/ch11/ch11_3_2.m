% % % % IRINA GRIGORESCU
% % % % DATE: 14-Nov-2016
% % % % 
% % % % CHAPTER 11.3.2
% % % % 
% % % % Lorentzian form - Problem 11.10
% % % % 

clear all; close all; clc
addpath ../../helpers

%% 
% Plot the sine and cosine functions
NSIZE = 1000; % points
xLim  =   20; % domain limit x \in(-xLim, xLim)
x     = linspace(-xLim, xLim, NSIZE); % domain of function
dx    = abs(x(2)-x(3));
dk    = 1/(NSIZE*dx);
k     = linspace(-(NSIZE*dk/2), (NSIZE*dk/2), NSIZE); % k-space domain
k(end)= 0;
k     = sort(k);

% k     = linspace(-(NSIZE*dk/2), 0, NSIZE/2); % k-space domain
% k     = [k linspace(abs(dk), NSIZE*dk/2, NSIZE/2)]; % k-space domain

%%
% Create function H(f) and its derivative
% h(t) = e^{-|t|/T2} and H(f) = (2a)/(a^2 + 4 pi^2 f^2), where a=1/T2

T2 = [50 80 100];

Hf = zeros(numel(T2), NSIZE);

for i = 1:numel(T2)
    a = 1/T2(i);
    Hf(i,:) = (2*a)./((a^2) + (4*pi^2).*(k.^2));
end

%%
% Plot Hf

figure,
for i = 1:numel(T2)
    plot(k(1:NSIZE), Hf(i,1:NSIZE))
    hold on
%     plot(k(NSIZE/2+1:end), Hf(i,NSIZE/2+1:end))
%     hold on
end
xlim([-0.5 0.5])
legend(strread(num2str(T2),'%s'));

%% 
% Create the function representing the area under H(f)

T2 = [50 80 100];
Hfint = zeros(numel(T2), NSIZE);

for i = 1:numel(T2)
    a = 1/T2(i);
    Hfint(i,:) = (sqrt(a)/pi) .* atan((2*pi/sqrt(a)).*k);
end

%%
% Plot Hfint

figure,
for i = 1:numel(T2)
    plot(k, Hfint(i,:))
    hold on
end
xlim([-2 2])
legend(strread(num2str(T2),'%s'));









