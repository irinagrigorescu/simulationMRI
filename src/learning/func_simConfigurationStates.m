function func_simConfigurationStates(N, T1, T2, ...
                   videoFlag)
% % % % IRINA GRIGORESCU
% % % % DATE CREATED: 11-07-2017
% % % % DATE UPDATED: 11-07-2017
% % % % 
% % % % Simulating the configuration states from
% % % % Extended Phase Graphs: 
% % % % Dephasing, RF Pulses, and Echoes - Pure and Simple
% % % % by Matthias Weigel
% % % % 
% % % % INPUT:
% % % %     N   = number of spins
% % % %     T1  = T1 relaxation time (ms)
% % % %     T2  = T2 relaxation time (ms)
% % % % 
% % % % OUTPUT: 
% % % %     - video
% % % % 

%% Make video if flag set to 1
if videoFlag.flag == 1
    % Make a frame by frame movie
    % need to set this up
    videoClass = VideoWriter([videoFlag.nameVideo, '.mp4'], ...
                             'MPEG-4'); % #need this for video

    frameRate = 25;                        % optional
    videoClass.set('FrameRate',frameRate); % optional
    open(videoClass)                       % #need this for video
end

%%
% % % % PREREQUISITES:
% We also know gamma = 4258 Hz/G (Hz/Gauss)
% 1 Tesla = 10^4  Gauss 
% 1 Gauss = 10^-4 Tesla = 100uT
gammabar =  42.58 * 1e+03; % 42.58 MHz/T = 42.58 * 10^3 (kHz / T )
grad     =  [0 ; 0 ; 50.00 * 1e-07]; %  0.05 G/mm  =  0.05 * 10^-4 ( T  / mm)
% % Positions of Spins in the z direction:
positionSpins = [zeros(1,N) ; zeros(1,N) ; linspace(-1, 1, N)]; % mm
% % Frequency induced due to gradient:
omegaGrad = (gammabar * 2*pi) .* (grad.' * positionSpins); % gamma*G*x (rad kHz)
% % The spins
M = repmat([1 0 0]', 1, N); % Lying in the transverse plane


% % % % FIGURE and its global details
figHand = figure('Position', [10,10,1200,800]); % #need this for video
dt = 0.1; % ms
Ntime = 50; % number of time points
c = linspace(1,2,N); % Colors for the magnetic moments tips

% % % % The event matrix
subplot(4,2,[7,8])
area(dt:dt:(Ntime*dt), ones(1,Ntime)); hold on
eventTip = stem(0, 1.2, 'k', 'filled');
set(gca, 'Ytick', [0 1]); 
set(gca, 'YtickLabel', {'0', 'G'}); 
set(gca, 'Xtick', []); 
xlim([0, (Ntime+1)*dt]); ylim([0, 2]);
xlabel('time'); ylabel('gradient'); 

% % % % The 3D plot
subplot(4,2,[1,6])
axisl = 1.1;
createAxis(axisl);
gline = plot3([0, 0], [0, 0], [-1, 1], 'k');
gline.LineWidth = 2;
% scatter3(0,0,1, 'kd', 'filled');
axis square; 
xlim([-1 1]); ylim([-1 1]); 
zlim([-axisl axisl]); 
view(35, 15); %view(130,15); %view(110, 10); %(110,10)

% % Plot the initial positions
MLinesStart = positionSpins; % start coordinate
MLinesEnd   = M;             % end coordinate
MLinesEnd(3, :) = positionSpins(3, :); 

vecLine3D = plot3([MLinesStart(1,:) ; MLinesEnd(1,:)], ...
                  [MLinesStart(2,:) ; MLinesEnd(2,:)], ...
                  [MLinesStart(3,:) ; MLinesEnd(3,:)], 'k');
% vecLine3D.Color(4) = 0.1;

hold on
vecTip3D  = scatter3(M(1, :), ...
                     M(2, :), ...
                     M(3, :) + positionSpins(3,:), ...
                     [], c, 'filled');

% % % Pause so I can start the simulation with a key
pause

% Plot for each dephasing
for i = 1:Ntime
    % % % % Delete previous plots
    delete(vecTip3D);
    delete(vecLine3D);
    delete(eventTip);
    
    % % % % Calculate dephasing due to gradient:
    phiGrad = omegaGrad * dt;
    
    % % % % Apply rotation due to gradient:
    % % % % for each spin
    for j = 1:N
        M(:, j) = Drel(dt, T1, T2)  * ...
                  Rotz(phiGrad(j)) * M(:, j) + ... % anti-clockwise rot
                  Drelz(dt, T1, 1);                % against all my beliefs 
    end                                            % :(
    
    % % % % Plot after dephasing
    MLinesStart = positionSpins; % start coordinate
    MLinesEnd   = M;             % end coordinate
    MLinesEnd(3, :) = positionSpins(3, :); 

    % % % % The event matrix
    subplot(4,2,[7,8])
    area(dt:dt:(Ntime*dt), ones(1,Ntime)); hold on
    eventTip = stem(dt*i, 1.2, 'k', 'filled');
    
    % % % % The 3D plot
    subplot(4,2,[1,6])
    vecLine3D = plot3([MLinesStart(1,:) ; MLinesEnd(1,:)], ...
                      [MLinesStart(2,:) ; MLinesEnd(2,:)], ...
                      [MLinesStart(3,:) ; MLinesEnd(3,:)], 'k');
	%vecLine3D.Color(4) = 0.1;
    
    hold on
    vecTip3D  = scatter3(M(1, :), ...
                         M(2, :), ...
                         M(3, :) + positionSpins(3,:), ...
                         [], c, 'filled');
    
    % % Save into video if flag set to 1
    if videoFlag.flag == 1
        frame = getframe(figHand);                 % #need this for video
        writeVideo(videoClass,frame)               % #need this for video
    end
    
    % % % % Pause for display purposes
    pause(0.1);
end

% % Close video handle if video flag set to 1
if videoFlag.flag == 1
    close(videoClass)                              % #need this for video
end


