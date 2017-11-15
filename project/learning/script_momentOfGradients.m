% % % % IRINA GRIGORESCU
% % % % Date created: 11-08-2017
% % % % 
% % % % This script is concerned with understanding the moment of gradients
% % % % It creates a gradient waveform and calculates it's first 3 moments.
% % % % 
close all
N = 100;
t = linspace(0,2,N)';
G1 = [ zeros(     N/5   ,1) ; ... % a tenth is              0
        ones((N-2*N/5)/3,1) ; ... % a third of the rest is  1
    -2.*ones((N-2*N/5)/3,1) ; ... % a third of the rest is -2
        ones((N-2*N/5)/3,1) ; ... % a third of the rest is  1
       zeros(     N/5   ,1)]; ... % a tenth is              0

G2 = [ zeros(     N/5   ,1) ; ... % a tenth is              0
        ones((N-2*N/5)/3,1) ; ... % a third of the rest is  1
    -2.*ones((N-2*N/5)/3,1) ; ... % a third of the rest is -2
     3.*ones((N-2*N/5)/3,1) ; ... % a third of the rest is  3
       zeros(     N/5   ,1)]; ... % a tenth is              0    
   
% Moments:
m0G1 = cumtrapz(t, G1);           % zeroth
m1G1 = cumtrapz(t, G1 .*  t);     % first
m2G1 = cumtrapz(t, G1 .* (t.^2)); % second

m0G2 = cumtrapz(t, G2);           % zeroth
m1G2 = cumtrapz(t, G2 .*  t);     % first
m2G2 = cumtrapz(t, G2 .* (t.^2)); % second

% % % % Plot the gradients and their moments
figure('Position', [100,400,1100,300]);
subplot(1,2,1)
plot(t, G1, 'k');      % G(t) 
hold on
plot(t, m0G1, 'r.');   % m0
plot(t, m1G1, 'g--');  % m1
plot(t, m2G1, 'b-');   % m2
grid on
ylim([-3.1 3.1])
xlabel('t(s)')
ylabel('G(t) and moments')
legend('G(t)', 'm_0', 'm_1', 'm_2')

subplot(1,2,2)
plot(t, G2, 'k');      % G(t) 
hold on
plot(t, m0G2, 'r.');   % m0
plot(t, m1G2, 'g--');  % m1
plot(t, m2G2, 'b-');   % m2
grid on
ylim([-3.1 3.1])
xlabel('t(s)')
ylabel('G(t) and moments')
legend('G(t)', 'm_0', 'm_1', 'm_2')