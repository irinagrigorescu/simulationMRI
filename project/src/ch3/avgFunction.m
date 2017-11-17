function avgFunction(n, m, omega, option)
% 
% Function that prints the cummulative average over the domain of a
% sinusoid or cosinusoid function 
% See Problem 3.2 from "Brown, MRI"
% 
% INPUT:
%   n = positive integer
%   m = positive integer
%   omega = frequency
%   option = a/b for the options of the problem
%            a  -> T  = 2*m*pi/omega
%            b  -> T >> 2*  pi/omega
% 
% OUTPUT:
%   Plots the cummulative integral over the domain
%   and the integral over the entire domain
% 
% Author: Irina Grigorescu, irina.grigorescu.15@ucl.ac.uk

% % % Prerequisites:
if nargin < 4
    option = 'a';
end

% % % % Parameters:
switch option
    case 'a'
        disp('Case a) chosen');
        T = 2*m*pi/omega;
    case 'b'
        disp('Case b) chosen');
        T = 2*pi/omega * 10000; 
end
N = 1000;
t = linspace(0, T, N);

% % % % Functions:
fsin    = sin(n*omega.*t);
fsinavg = cumtrapz(t,fsin).*(1./t);
fcos    = cos(n*omega.*t);
fcosavg = cumtrapz(t,fcos).*(1./t);

% % % % Plotting:
figure('Position',[10,10,1200,700])
% 1. sin
subplot(2,2,1)
plot(t, fsin)
xlim([min(t) max(t)])
xlabel('t'); ylabel('sin(n \omega t)'); 
title('f(t) = sin(n \omega t)');

% 2. Average
subplot(2,2,2)
plot(t, fsinavg)
xlim([min(t) max(t)])
xlabel('t'); ylabel('1/t \int_0^t dt sin(n \omega t)'); 
title(['<f(t)> \equiv 1/T \int_0^T dt f(t) = ', ...
       '1/T \int_0^T dt sin(n \omega t) = ', num2str(fsinavg(end))]);

% 3. cos
subplot(2,2,3)
plot(t, fcos)
xlim([min(t) max(t)])
xlabel('t'); ylabel('cos(n \omega t)'); 
title('f(t) = cos(n \omega t)');

% 4. Average
subplot(2,2,4)
plot(t, fcosavg)
xlim([min(t) max(t)])
xlabel('t'); ylabel('1/t \int_0^t dt cos(n \omega t)'); 
title(['<f(t)> \equiv 1/T \int_0^T dt f(t) = ', ...
       '1/T \int_0^T dt cos(n \omega t) = ', num2str(fcosavg(end))]);




