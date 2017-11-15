% IRINA GRIGORESCU
% 
% The Spiral form of EPI
% Playing around with it
% 

clear all; close all

% Prerequisites
gamma    = 2.68*1E08;      %(rad/s/T)
gammabar = gamma./(2*pi);  %(1/s/T)

% Imaging requirements
FOV = 0.3;                    % 300mm = 0.3m
Nx = 192; Ny = 192;           % Nx = Ny = 192
dx = FOV/Nx;
Ts = 0.06;                    % 6ms = 0.06s
fBW = 50000;                  % BW frequency = 50kHz = 50,000 Hz
dt = 1/fBW;                   % dwell time dt = 1/BW
Ns = Ts/dt;                   % number of samples
Nr = 48; Nt = ceil(Ns/Nr);
Gmax = 1 / (gammabar*dt*FOV); 

% K-space trajectory
Nshot = 1;
t = 0:dt:Ts-dt;                            % Ns time points from 0 to 6 ms
alpha1 = pi/(gamma*Nt*Nr*dt*dx);
alpha2 = 2*pi/(Nt*dt);                  % thetaMax = pi N / Nshot
                                           % thetaMax = alpha*Ts
                                           % alpha = (pi N) / (Ts Nshot)
kx  = gamma * alpha1 .* (t.^2.1) .* cos(alpha2 .* t);
ky  = gamma * alpha1 .* (t.^2.1) .* sin(alpha2 .* t);

% Gradient trajectory
Gx = alpha1 .* cos(alpha2.*t) - alpha1 .* alpha2 .* t .* sin(alpha2.*t);
Gy = alpha1 .* sin(alpha2.*t) + alpha1 .* alpha2 .* t .* cos(alpha2.*t);

% Plot
figure; hold on
subplot(2,2,1)
scatter(kx, ky, '.'); 
xlabel('k_x (m^{-1})'); ylabel('k_y (m^{-1})'); 
axis equal

subplot(2,2,2)
scatter(t, Gx.*1E+03, '.'); 
xlabel('s'); ylabel('G_x (mT)'); 

subplot(2,2,4)
scatter(t, Gy.*1E+03, '.'); 
xlabel('s'); ylabel('G_y (mT)'); 










%%
% Functions
kx_traj = @(a1, a2, t) gammabar .* a1 .* t .* sin(a2.*t);
ky_traj = @(a1, a2, t) gammabar .* a1 .* t .* cos(a2.*t);
gx_traj = @(a1, a2, t) a1 .* sin(a2.*t) + a1 * a2 .* t .* cos(a2 .* t);
gy_traj = @(a1, a2, t) a1 .* sin(a2.*t) - a1 * a2 .* t .* sin(a2 .* t);

% Parameters
N = 200;
alpha1 = deg2rad( 1);
alpha2 = deg2rad(35);
t = linspace(0,50,N);

% Calculate trajectory
kxtr = kx_traj(alpha1, alpha2, t);
kytr = ky_traj(alpha1, alpha2, t);
gxtr = gx_traj(alpha1, alpha2, t);
gytr = gy_traj(alpha1, alpha2, t);

maxVK = max( max(kxtr), max(kytr) );
minVK = min( min(kxtr), min(kytr) );

maxVG = max( max(gxtr), max(gytr) );
minVG = min( min(gxtr), min(gytr) );


% Plot trajectory
figure(1)%('Position', [10,10, 1200, 800])
c = 'b';
nr = 4; nc = 5;
for i = 1:N
    
    % Gradient waveforms - Gx
    subplot(nr, nc, [1, 8])
    plot(t(1:i), gxtr(1:i), [c, '.-']); hold on
    grid on
    xlim([0, max(t)]);   
    ylim([minVG, maxVG]);
    xlabel('time'); ylabel('G_x');
    title('Gradient waveform (G_x)')
    
    % Gradient waveforms - Gy
    subplot(nr, nc, [11, 18])
    plot(t(1:i), gytr(1:i), [c, '.-']); hold on
    grid on
    xlim([0, max(t)]);   
    ylim([minVG, maxVG]);
    xlabel('time'); ylabel('G_y');
    title('Gradient waveform (G_y)')
    
    % Gradient vector - G
    subplot(nr, nc, [14, 20])
    vec1 = plot([0, gxtr(i)], [0, gytr(i)], [c, '.-']); hold on
    vec2 = plot([0, gxtr(i)], [0, 0], 'r'); 
    vec3 = plot([0, 0], [0, gytr(i)], 'g'); 
    grid on   
    xlim([minVG, maxVG]);
    ylim([minVG, maxVG]);
    xlabel('G_x'); ylabel('G_y');
    axis square
    title({'Gradient trajectory'})
    
    % K-space trajectory
    subplot(nr, nc, [4, 10])
    plot(kxtr(1:i), kytr(1:i), [c, '.-']); hold on
    xlim([minVK, maxVK]);   
    ylim([minVK, maxVK]);
    grid on
    xlabel('k_x'); ylabel('k_y');
    axis square
    title({'K-space trajectory for', ...
           ['\alpha_1 = ', num2str(rad2deg(alpha1)), ...
            '^o and \alpha_2 = ', num2str(rad2deg(alpha2)), '^o']});
    
    drawnow; %pause(0.01)
    
    if (i ~= N)
        delete(vec1);
        delete(vec2);
        delete(vec3);
    end
    
end
    
    