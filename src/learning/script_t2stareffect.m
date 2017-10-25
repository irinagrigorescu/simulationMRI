% % Irina Grigorescu
% This script looks at the T2star effect
% by simulating the signals coming from magnetic moments of different
% off-resonance frequencies and then the net signal coming from the sum
% over them.

% Note: off-resonance values reflect 1/T2'
% By choosing T2star and T2, you can calculate T2'

% Number of timepoints
N  = 1000;
% Relaxation rate
T2  =  30; % ms
T2s =  15; % ms
T2p =  T2*T2s/(T2-T2s); % ms
R2p = 1/T2p;
w0  =   50; % kHz

% Signal function
sig = @(t,dw) exp(-t./T2) .* cos((w0+dw).*t);

% A few off-resonance values
dws = -2*R2p + 4*R2p.*randn(1,1000);
idx = [1:5 20:25 1000];
% Time domain
t  = linspace(0, 150, N);

% Calculate signal
Msig = zeros(size(dws,2), N);

% Calculate signal for each off-resonance value
for i = 1:size(dws,2)
    Msig(i,:) = sig(t,dws(i));
end

% % % % PLOT
figure

% Plot the envelope
plot(t, exp(-t./T2), 'g', 'LineWidth', 1); hold on
plot(t, exp(-t./T2s), 'r', 'LineWidth', 1); 
ylim([-1 1]);
xlabel('time (ms)'); ylabel('signal (a.u.)');
legend('exp(-t/T_2)', 'exp(-t/T_2^*)');
pause

% Plotting some of the individual spins
for i = 1:size(idx,2)
    tp = plot(t, Msig(idx(i),:)); 
    pause(0.1)
    delete(tp);
end

% Plot their net sum
plot(t, 1/size(dws,2).*sum(Msig,1), 'k', 'LineWidth', 1)




