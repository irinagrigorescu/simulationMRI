% % % % 
% % % % IRINA GRIGORESCU
% % % % Date created: 29-06-2017
% % % % Date updated: 29-06-2017
% % % % 
% % % % This script looks at selective excitation.
% % % % Only spins with a certain range of frequencies will be excited.
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % Website: http://mrsrl.stanford.edu/~brian/bloch/
% % % % F-1. Basic Selective Excitation
% % % % 

% % % Prerequisites:
org = [0 0 0]; % origin point
xax = [1 0 0]; yax = [0 1 0]; zax = [0 0 1];
xyzax = [xax ; yax ; zax];

% % 0. Tissue with:
T1 = 600;
T2 = 100;

% % % % % % % % % % % No gradient, just RF pulses
% % 1) Two consecutive RF pulses of 45deg, 2.3ms apart
% RF pulse:
RF       = struct;
RF.alpha = pi/4;     %  flip angle
RF.phi   = 0;        % phase angle
% Time delay:
tau = 2.3;         % in ms

% % 1.2) Plot the transverse magnetisation (magnitude and phase) immediately
% after the second RF pulse, over the resonant frequency range 
% [-500Hz, 500Hz]. 
% The excitation is selective in resonant frequency
freq = -0.5 : 0.005 : 0.5; % Hz frequencies
M0   = [0 0 1]'; % relaxed magnetisation vector 

Mxyz = zeros(3,length(freq));

for fr = 1:length(freq)
    % Apply first pulse
    M = Rotz(-RF.phi) * Rotx(-RF.alpha) * Rotz(RF.phi) * M0;
    
    % After first pulse free precession happens
    betaFreePrecess = 2*pi .* freq(fr) .* tau; % dw*t=2pi*df*t
    M = Rotz(-betaFreePrecess) * Drel(tau, T1, T2) * M + ...
        Drelz(tau, T1, 1);
    
    % Apply next RF pulse
    M = Rotz(-RF.phi) * Rotx(-RF.alpha) * Rotz(RF.phi) * M;
    
    Mxyz(1, fr) = M(1);
    Mxyz(2, fr) = M(2);
    Mxyz(3, fr) = M(3);
end

% Plot for each frequency the signal magnitude and phase
figure
subplot(2,1,1)
plot(freq.*1000, abs(Mxyz(1, :) + 1i.*Mxyz(2, :)))
xlabel('Frequency (Hz)');
ylabel('Signal (fraction of M_0)');
grid on;
title('A 2 RF pulse experiment')

subplot(2,1,2)
plot(freq.*1000, atan2(Mxyz(2, :), Mxyz(1, :)))
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
grid on;

% % % % % % % % % % % 
% % % % % % % % % % % 
% % % % % % % % % % % We turn on a gradient of strength 0.1G/cm in the x
% direction for the whole excitation. 
% We also know gamma = 4258 Hz/G (Hz/Gauss)
% 1 Tesla = 10^4  Gauss 
% 1 Gauss = 10^-4 Tesla = 100uT
xpos   = -20 : 0.01 : 20; % mm in position
gamma = 4.258;          % kHz/G
grad  = 0.005;          % 0.1 G/cm = 0.01G/mm
M0   = [0 0 1]'; % relaxed magnetisation vector 
tau = 2.3;

Mxyz = zeros(3, length(xpos));
for x = 1:length(xpos)
    % Apply first pulse
    M = Rotz(-RF.phi) * Rotx(-RF.alpha) * Rotz(RF.phi) * M0;
    
    % After first pulse free precession happens (no free precession)
    betaFreePrecess = 2*pi .* 0.2 .* tau; % dw*t=2pi*df*t
    M = Rotz(-betaFreePrecess) * Drel(tau, T1, T2) * M + ...
        Drelz(tau, T1, 1);
    
    % Rotation due to gradient
    phiGrad = (gamma * 2*pi) * grad * xpos(x) * tau; % gamma*G*x*t
    M = Rotz(-phiGrad) * M;
    
    % Apply next RF pulse
    M = Rotz(-RF.phi) * Rotx(-RF.alpha) * Rotz(RF.phi) * M;
    
    % Apply refocusing gradient
    phiGrad = (gamma * 2*pi) * -(0.5*grad) * xpos(x) * tau; % gamma*G*x*t
    M = Rotz(-phiGrad) * M;
    
    Mxyz(1, x) = M(1);
    Mxyz(2, x) = M(2);
    Mxyz(3, x) = M(3);
end

% Plot for each frequency the signal magnitude and phase
figure
subplot(2,1,1)
plot(xpos, abs(Mxyz(1, :) + 1i.*Mxyz(2, :)))
xlabel('Position (mm)');
ylabel('Signal (fraction of M_0)');
grid on;
title({'With a gradient on between the 2 RF pulses', ...
       'and a refocusing gradient of half area'});

subplot(2,1,2)
plot(xpos, atan2(Mxyz(2, :), Mxyz(1, :)))
xlabel('Position (mm)');
ylabel('Phase (rad)');
grid on;




% % %%
% % % Plot before and after pulse:
% % figure
% % 
% % % Plot xyz
% % for i = 1:3
% %     plot3([-xyzax(i,1), xyzax(i,1)], ...
% %           [-xyzax(i,2), xyzax(i,2)], ...
% %           [-xyzax(i,3), xyzax(i,3)], 'k--'), hold on
% % end
% % % Before pulse:
% % plot3([org(1), M0(1)], [org(2), M0(2)], [org(3), M0(3)], 'b'), hold on
% % scatter3(M0(1), M0(2), M0(3), 'bd'), hold on
% % xlabel('x'); ylabel('y'); zlabel('z'); 
% % grid on
% % axis equal
% % % B1 field:
% % B1 = [1 0 0]';             
% % B1 = Rotz(-RF.phi) * B1;
% % plot3([org(1),  B1(1)], [org(2),  B1(2)], [org(3), B1(3)], 'gd--');
% % view(0,90)
% % % Apply first pulse:
% % for i = 1:10
% %     M = Rotz(-RF.phi) * Rotx(-(RF.alpha).*i/10) * Rotz(RF.phi) * M0;
% %     vec = plot3([org(1),  M(1)], [org(2),  M(2)], [org(3), M(3)], 'r');
% %     tip = scatter3(M(1), M(2), M(3), 'rd');
% %     pause(1)
% %     if i ~= 10
% %         delete(vec)
% %     end
% % end
% % M0 = M;
% % % Free precession:
% % for i = 1:10
% %     betaFreePrecess = 2*pi .* freq(102) .* tau .* i / 10; % dw*t=2pi*df*t
% %     
% %     M = Rotz(-betaFreePrecess) * Drel(tau.*i/10, T1, T2) * M0 + ...
% %         Drelz(tau.*i/10, T1, 1);
% %     
% %     vec = plot3([org(1),  M(1)], [org(2),  M(2)], [org(3), M(3)], 'r');
% %     tip = scatter3(M(1), M(2), M(3), 'rd');
% %     pause(1)
% %     if i ~= 10
% %         delete(vec)
% %     end
% % end
% % M0 = M;
% % % Apply second pulse:
% % for i = 1:10
% %     M = Rotz(-RF.phi) * Rotx(-(RF.alpha).*i/10) * Rotz(RF.phi) * M0;
% %     vec = plot3([org(1),  M(1)], [org(2),  M(2)], [org(3), M(3)], 'r');
% %     tip = scatter3(M(1), M(2), M(3), 'rd');
% %     pause(1)
% %     if i ~= 10
% %         delete(vec)
% %     end
% %     if i == 10
% %         plot3([org(1),  M(1)], [org(2),  M(2)], [org(3), M(3)], 'k');
% %         scatter3(M(1), M(2), M(3), '*');
% %     end
% % end

