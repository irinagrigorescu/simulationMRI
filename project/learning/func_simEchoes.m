function func_simEchoes(N, T1, T2, dws, ...
                   rfa, pha, tau, Ttotal, ...
                   videoFlag)
% % % % IRINA GRIGORESCU
% % % % DATE CREATED: 28-06-2017
% % % % DATE UPDATED: 28-06-2017
% % % % 
% % % % Simulating echoes in a graphical way
% % % % 
% % % % INPUT:
% % % %     N   = number of spins
% % % %     T1  = T1 relaxation time (ms)
% % % %     T2  = T2 relaxation time (ms)
% % % %     dws = off-resonance frequencies for the N spins (kHz)
% % % %     rfa = flip  angles (deg)
% % % %     pha = phase angles (deg)
% % % %     tau = time between rf pulses (ms)
% % % %     Ttotal = total time of simulation (ms)
% % % %     videoFlag = (struct) with:
% % % %                .flag      = 1 if I want a video
% % % %                .nameVideo = name of the video (string)
% % % % 
% % % % OUTPUT: 
% % % %     - video
% % % % 

isRFevents   = 1;
is3DView     = 1;
is2DView     = 1;
isSignal     = 1;
isComponents = 0;

%% Make video if flag set to 1
if videoFlag.flag == 1
    % Make a frame by frame movie
    % need to set this up
    videoClass = VideoWriter([videoFlag.nameVideo, '.mp4'], ...
                             'MPEG-4'); % #need this for video

    frameRate = 2;                        % optional
    videoClass.set('FrameRate',frameRate); % optional
    open(videoClass)                       % #need this for video
end

%% Create event matrix for plotting

% Prerequisites:
Npulses = length(rfa); % how many pulses we have
dt  = 0.04;             % how finely grained you want your simulation

% See which indices correspond to negative/positive off-resonance freq
indicesN = find(dws<=0);
indicesP = find(dws >0);

% Events matrix with:
Tpoints = ceil(Ttotal/dt); % number of time points in simulation
% % % Events: 1 = timepoint | 2 = RF.FA | 3 = RF.PhA
RFevents = zeros(3, Tpoints);
% Populate event matrix with the timpoints (1st line)
RFevents(1, :) = linspace(0, Ttotal, Tpoints); 
% Populate eventmatrix with the rf pulses and phase angles
for i = 1:Npulses
	RFevents(2, floor(sum(tau(1:i))/dt)+1) =  rfa(i); % 2nd line
    RFevents(3, floor(sum(tau(1:i))/dt)+1) =  pha(i); % 3rd line
end

% Keep position of all of them in M (3xN)
M   = repmat([0 0 1]', 1, N);    % The m vectors
% The 0,0,0 -> Mxyz lines
MLines = zeros(3, 2*N);          % The vector lines (start from 0,0,0)
MLines(:, 2:2:end) = M(:, :);    % end coordinate   (end at M_x,y,z)

% % % % Figure and its global details:
figHand = figure('Position', [10,10,1400,1200]); % #need this for video
nr = 5; nc = 11;         % subplots
eventsPos = [ 1,  5];    % event matrix
signalPos = [ 7, 11];    % signal
simulMxyz = [18, 22];    % The x,y,z components
simul3DPos = [  12, 49]; % the 3D view 
simul2DPos = [18,55]; %[  29, 55]; % the 2D projection view
c = linspace(1,2,N);     % Colors for the magnetic moments tips
axisl = 1.1;             % axis limits in the 3D/2D plots

% % %
% 1st subplot is the Events
if isRFevents
    subplot(nr,nc,eventsPos)
    % The phase angles
    scatter(RFevents(1,:), RFevents(3,:), 'ro', 'filled'), hold on
    % The flip angles
    stem(   RFevents(1,:), RFevents(2,:), 'b')
    % Details about the plot
    grid on
    xlim([-0.1, Ttotal])
    xlabel('t (ms)'); ylabel('RF angle (deg)')
    title('RF Events');
end

% % %
% 2nd subplot is the Magnetic moments (the 3D plot)
if is3DView
    subplot(nr,nc,simul3DPos)
    % Create the coordinate system axis
    createAxis(axisl); 
    % The vector tip
    vecTip3D1 = scatter3(M(1, indicesN), ...
                         M(2, indicesN), ...
                         M(3, indicesN), ...
                         [], 'ob', 'filled'); hold on
    vecTip3D2 = scatter3(M(1, indicesP), ...
                         M(2, indicesP), ...
                         M(3, indicesP), ...
                         [], 'or', 'filled'); hold on
    % The vector line
    vecLine3D = plot3(MLines(1, :), ...
                      MLines(2, :), ...
                      MLines(3, :), 'k' );
    % Plot details
    axis square; %xlim([-1 1]); ylim([-1 1]); zlim([-1 1]); 
    xlabel('M_x') ; ylabel('M_y'); zlabel('M_z');
    view(110, 10); %(110,10)
    title('3D View');
end

% % %
% 3rd subplot is the Magnetic moments (the 2D plot)
if is2DView
    subplot(nr,nc,simul2DPos)
    % The coordinate system lines
    plot([0 0], [1 -1], 'k--'), hold on, plot([1 -1], [0 0], 'k--');
    % The vector tip
    vecTip2D1 = scatter(M(1, :), M(2, :), [], 'ob', 'filled'); hold on
    vecTip2D2 = scatter(M(1, :), M(2, :), [], 'or', 'filled'); hold on
    % The vector line
    vecLine2D = plot(MLines(1, :), MLines(2, :), 'k' );
    vecLineMNet = plot([0 M(1,1)], [0 M(2,1)], 'c', 'LineWidth', 2);
    % Plot details
    grid on
    axis equal
    xlim([-axisl axisl]); ylim([-axisl axisl])
    xlabel('M_x') ; ylabel('M_y'); 
    title('View of Transverse Plane');
end

% % %
% 4th subplot is the Signal
if isSignal
    subplot(nr,nc,signalPos)
    % Add T2 exponential decay curve
    plot(RFevents(1,:), exp(-RFevents(1,:)/T2)), hold on
    % Add the RF pulses
    stem(RFevents(1,:), RFevents(2,:)~=0, 'k--', 'marker', 'none')
    % Plot details
    grid on
    xlim([-0.1, Ttotal]); ylim([0, 1]);
    ylabel('|M_{xy}|'); xlabel('t (ms)');
    legend('-e^{-t/T_2}');
    title('Signal Evolution');
end

% % %
% 5th subplot is the components
if isComponents
    subplot(nr,nc,simulMxyz)
    netM = sum(M,2)./N;
    scatter(RFevents(1,i),netM(1),'b.'); hold on
    scatter(RFevents(1,i),netM(2),'r.');
    scatter(RFevents(1,i),netM(3),'g.');
    grid on
    xlim([0, Ttotal]); ylim([-1, 1]);
    ylabel('M_i'); xlabel('t (ms)');
    legend('M_x','M_y','M_z')
    title('Magnetic moment vector components');
end

% % Save into video if flag set to 1
if videoFlag.flag == 1
    frame = getframe(figHand);                 % #need this for video
    writeVideo(videoClass,frame)               % #need this for video
end

% Pause before simulation
pause

% Do accelerated dephasing after 2nd rf pulse
isSecondRFPulse = -1; % added now irina
Mnet = zeros(3,1);
liniarIncreaseDW = 0; %0.001;

% Start simulation
for i = 1:Tpoints
    
    figure(1)
    
    % Delete history of vectors
    if is3DView
        delete(vecTip3D1); delete(vecTip3D2); 
        delete(vecLine3D);
    end
    if is2DView
        delete(vecTip2D1); delete(vecTip2D2);
        delete(vecLine2D); delete(vecLineMNet);
    end
    
    % SUBPLOT 1 - Plot events
    if isRFevents
        subplot(nr,nc,eventsPos)
        event = stem(RFevents(1,i), 90, 'k');
    end
    
    % DO RF PULSE
    if RFevents(2,i) ~= 0
        isSecondRFPulse = isSecondRFPulse + 1; % added now irina
        
        for j = 1:10                        
            % Update the tips during RF Pulse
            Mtemp = Rotz( deg2rad(RFevents(3,i))) * ...
                     Rotx(-deg2rad(RFevents(2,i)*j/10)) * ...
                      Rotz(-deg2rad(RFevents(3,i))) * M;
            % Update the lines during RF Pulse (they always start from 0,0,0)
            MLinesTemp = MLines;
            MLinesTemp(:, 2:2:end) = Mtemp(:, :);  % end   coordinate
                  
            % % % THE 3D plot
            if is3DView
                subplot(nr,nc,simul3DPos)
                % Plot the lines
                vecLine3D = plot3(MLinesTemp(1, :), ...
                                  MLinesTemp(2, :), ...
                                  MLinesTemp(3, :), 'k' );
                vecLine3D.Color(4) = 0.5;
                % Plot the tips
                vecTip3D1  = scatter3(Mtemp(1, indicesN), ...
                                     Mtemp(2, indicesN), ...
                                     Mtemp(3, indicesN), ...
                                     [], 'ob', 'filled');
                vecTip3D2  = scatter3(Mtemp(1, indicesP), ...
                                      Mtemp(2, indicesP), ...
                                      Mtemp(3, indicesP), ...
                                      [], 'or', 'filled');
                hold on
            end
                       
            % % % THE 2D plot
            if is2DView
                subplot(nr,nc,simul2DPos)
                % Plot the lines
                vecLine2D = plot(MLinesTemp(1, :), ...
                                 MLinesTemp(2, :), 'k' );
                vecLine2D.Color(4) = 0.5;
                % Plot the tips
                vecTip2D1 = scatter(Mtemp(1, indicesN), ...
                                    Mtemp(2, indicesN), ...
                                    [], 'ob', 'filled');
                vecTip2D2 = scatter(Mtemp(1, indicesP), ...
                                    Mtemp(2, indicesP), ...
                                    [], 'or', 'filled');
            end
                           
            % % % The components plot
            if isComponents
                subplot(nr,nc,simulMxyz)
                netM = sum(Mtemp,2)./N;
                scatter(RFevents(1,i),netM(1),'b.'); hold on
                scatter(RFevents(1,i),netM(2),'r.');
                scatter(RFevents(1,i),netM(3),'g.');
                grid on
                xlim([0, Ttotal]); ylim([-1, 1]);
            end

            % Pause
            pause(0.1)
            % % Save into video if flag set to 1
            if videoFlag.flag == 1
                frame = getframe(figHand);                 % #need this for video
                writeVideo(videoClass,frame)               % #need this for video
            end
            
            % Delets
            if is3DView
                delete(vecTip3D1); delete(vecTip3D2); 
                delete(vecLine3D);
            end
            if is2DView
                delete(vecTip2D1); delete(vecTip2D2);
                delete(vecLine2D); delete(vecLineMNet);
            end
        end
        % Update the tips
        M = Rotz( deg2rad(RFevents(3,i))) * ...
                Rotx(-deg2rad(RFevents(2,i))) * ...
                    Rotz(-deg2rad(RFevents(3,i))) * M ;
        % Update the lines (they always start from 0,0,0)
        MLines(:, 2:2:end) = M(:, :);  % end   coordinate
    else
	% DEPHASE if not RF pulse events
    
        % If it's the second rf pulse
        if isSecondRFPulse == 1
            dws = dws + liniarIncreaseDW;
        end
    
        % Start dephasing for all magnetic moment vectors
        for j = 1:N
            % Update the tips
            M(:, j) = Drel(dt, T1, T2) * ...
                      Rotz(-2*pi*dt*dws(j)) * M(:, j) + ...
                      Drelz(dt, T1, 1);
            % Update the lines (they always start from 0,0,0)
            MLines(:, 2:2:end) = M(:, :);  % end   coordinate
        end
    end
    
    % % % The 3D plot
    if is3DView
        subplot(nr,nc,simul3DPos)
        % Plot the lines
        vecLine3D = plot3(MLines(1, :), ...
                          MLines(2, :), ...
                          MLines(3, :), 'k' );
        vecLine3D.Color(4) = 0.5;
        % Plot the tips
        vecTip3D1  = scatter3(M(1, indicesN), ...
                              M(2, indicesN), ...
                              M(3, indicesN), ...
                              [], 'ob', 'filled'); 
        vecTip3D2  = scatter3(M(1, indicesP), ...
                              M(2, indicesP), ...
                              M(3, indicesP), ...
                              [], 'or', 'filled'); 
        hold on
    end

    % % % The 2D plot
    if is2DView
        subplot(nr,nc,simul2DPos)
        % Plot the lines
        vecLine2D = plot(MLines(1, :), ...
                         MLines(2, :), 'k' );
        vecLine2D.Color(4) = 0.5;
        vecLineMNet = plot([0 Mnet(1)], [0 Mnet(2)], 'c', 'LineWidth', 2);
        
        % Plot the tips
        vecTip2D1 = scatter(M(1, indicesN), M(2, indicesN), [], 'ob', 'filled');
        vecTip2D2 = scatter(M(1, indicesP), M(2, indicesP), [], 'or', 'filled');

        % % % Signal
        if isSignal
            subplot(nr,nc,signalPos)
            scatter(RFevents(1,i), sqrt((sum(M(1,:))./N).^2 + ...
                                        (sum(M(2,:))./N).^2), 'k.')
            hold on
            grid on
            xlim([0, Ttotal])
        end
    end
    
    % % % The components
    if isComponents
        subplot(nr,nc,simulMxyz)
        netM = sum(M,2)./N;
        scatter(RFevents(1,i),netM(1),'b.'); hold on
        scatter(RFevents(1,i),netM(2),'r.');
        scatter(RFevents(1,i),netM(3),'g.');
        grid on
        xlim([0, Ttotal]); ylim([-1, 1]);
        ylabel('M_i'); xlabel('t (ms)');
    end
    
    % Pause
    drawnow
    
    % % Save into video if flag set to 1
    if videoFlag.flag == 1
        frame = getframe(figHand);                 % #need this for video
        writeVideo(videoClass,frame)               % #need this for video
    end

    Mphase = rad2deg(angle(M(1,:)+1i.*M(2,:))+pi);
    Mmagn  = abs(M(1,:)+1i.*M(2,:));
    Mnet   = mean(M,2);
%     figure(2)
%     subplot(1,2,1)
%     hist(Mphase,N/4)
%     subplot(1,2,2)
%     hist(Mmagn,N/4)
    
    % Delete RF events
    if isRFevents
        figure(1)
        delete(event)
    end
    
end



% % Close video handle if video flag set to 1
if videoFlag.flag == 1
    close(videoClass)                              % #need this for video
end

% % % end % FUNCTION END

