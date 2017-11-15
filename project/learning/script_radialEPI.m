% IRINA GRIGORESCU
% 
% The Radial form of EPI
% Playing around with it
% 

clear all; close all

gammabar = 42.58; % MHz T^-1

% Functions
kx_traj = @(k, th) k .* cos(th);
ky_traj = @(k, th) k .* sin(th);
gx_traj = @(G, th) G .* cos(th);
gy_traj = @(G, th) G .* sin(th);

% Parameters
N =  20;
G =   2; % T/m
theta = deg2rad(0:25:179);
t  = linspace(0,50,N);
dt = t(2) - t(1);

maxVG = 1.5*G;             minVG = -maxVG;
maxVK = gammabar*G*max(t); minVK = -maxVK;

% Plot trajectory
figure('Position', [10,10, 1200, 800])
c = 'b';
nr = 4; nc = 5;

% Calculate gradient amplitude
gx = gx_traj(G, theta);
gy = gy_traj(G, theta);

for j = 1:size(theta,2)
    
    % Calculate trajectory
    kxtr = kx_traj(gammabar.*(-G).*t, theta(j));
    kytr = ky_traj(gammabar.*(-G).*t, theta(j));

    fprintf('j = %d => theta = %3.2f G_x = %3.2f G_y = %3.2f\n\n', ...
             j, rad2deg(theta(j)), gx(j), gy(j));

    for i = 1:N

        % Gradient waveforms - Gx
        subplot(nr, nc, [1, 8])
        if i == 1
            plot(t, repmat(gx(j), size(t,2),1), c); hold on
        end
        plot(t(i), gx(j), [c, 'o-']); 
        grid on
        xlim([min(t)-dt, max(t)+dt]);   
        ylim([minVG, maxVG]);
        xlabel('time'); ylabel('G_x');
        title('Gradient waveform ')

        % Gradient waveforms - Gy
        subplot(nr, nc, [11, 18])
        if i == 1
            plot(t, repmat(gy(j), size(t,2),1), c); hold on
        end
        plot(t(i), gy(j), [c, 'o-']); 
        grid on
        xlim([min(t)-dt, max(t)+dt]);   
        ylim([minVG, maxVG]);
        xlabel('time'); ylabel('G_y');

        % K-space trajectory
        subplot(nr, nc, [9, 15])
        plot(kxtr(1:i), kytr(1:i), [c, '.-']); hold on
        xlim([minVK, maxVK]);   
        ylim([minVK, maxVK]);
        grid on
        xlabel('k_x'); ylabel('k_y');
        axis square
        title({'K-space trajectory for', ...
               ['\theta = ', num2str(rad2deg(theta(j))), '^o']});

        pause(0.001)

    end
end
    