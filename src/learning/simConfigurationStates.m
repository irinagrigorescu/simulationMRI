% % function simConfigurationStates()
% % % % IRINA GRIGORESCU
% % % % DATE CREATED: 11-07-2017
% % % % DATE UPDATED: 11-07-2017
% % % % 
% % % % Simulating the configuration states from
% % % % Extended Phase Graphs: 
% % % % Dephasing, RF Pulses, and Echoes - Pure and Simple
% % % % by Matthias Weigel

% clear all
% close all
addpath ../helpers/

% % % % Prerequisites:
% We also know gamma = 4258 Hz/G (Hz/Gauss)
% 1 Tesla = 10^4  Gauss 
% 1 Gauss = 10^-4 Tesla = 100uT
gammabar =  42.58 * 1e+03; % 42.58 MHz/T = 42.58 * 10^3 (kHz / T )
grad     =  50.00 * 1e-07; %  0.05 G/mm  =  0.05 * 10^-4 ( T  / mm)

% % % % Choose Yourself:
% % Relaxation Terms:
T1 =  10000;
T2 =  10000;

% % Spins:
N = 1000;
positZ = linspace(0, 1, N); % mm

% % Frequency induced due to gradient:
omegaGrad = (gammabar * 2*pi) * grad * positZ; % gamma*G*x (rad kHz)

% % The spins
M = repmat([0 1 0]', 1, N); % Lying in the transverse plane


% % % % Plot
figure
axisl = 1.1;
createAxis(axisl);
axis square; zlim([0 axisl]); 
view(90,0); %view(110, 10); %(110,10)

% % Colors for the magnetic moments tips
c = linspace(1,2,N); 

% % Plot the initial positions
MLines = zeros(3, 2*N);
MLines(3, 1:2:end) = positZ(:); % start coordinate
MLines(:, 2:2:end) = M(:, :);   % end   coordinate
MLines(3, 2:2:end) = MLines(3, 2:2:end) + positZ; % to separate them

vecLine3D = plot3(MLines(1, :), ...
                  MLines(2, :), ...
                  MLines(3, :), 'k' );
hold on
vecTip3D  = scatter3(M(1, :), ...
                     M(2, :), ...
                     M(3, :) + positZ, ...
                     [], c, 'filled');

dt = 0.1; % ms
% % % Pause so I can start the simulation with a key
pause
% Plot for each dephasing
for i = 1:200
    % % % % Delete previous plots
    delete(vecTip3D);
    delete(vecLine3D);
    
    % % % % Calculate dephasing due to gradient:
    phiGrad = omegaGrad * dt;
    
    % % % % Apply rotation due to gradient:
    % % % % for each spin
    for j = 1:N
        M(:, j) = Drel(dt, T1, T2)  * ...
                  Rotz(-phiGrad(j)) * M(:, j) + ...
                  Drelz(dt, T1, 1);
    end
    
    % % % % Plot after dephasing
    MLines(3, 1:2:end) = positZ(:); % start coordinate
    MLines(:, 2:2:end) = M(:, :);   % end   coordinate
    MLines(3, 2:2:end) = MLines(3, 2:2:end) + positZ; % to separate them

    vecLine3D = plot3(MLines(1, :), ...
                      MLines(2, :), ...
                      MLines(3, :), 'k' );
    hold on
    vecTip3D  = scatter3(M(1, :), ...
                         M(2, :), ...
                         M(3, :) + positZ, ...
                         [], c, 'filled');
    
    % % % % Pause for display purposes
    pause(0.1);
end













