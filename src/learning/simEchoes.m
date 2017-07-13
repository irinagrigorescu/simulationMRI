% % function simEchoes()
% % % % IRINA GRIGORESCU
% % % % DATE CREATED: 28-06-2017
% % % % DATE UPDATED: 28-06-2017
% % % % 
% % % % Simulating echoes in a graphical way
% % % % 
% % % % 

% clear all
% close all
addpath ../helpers/

% % % % Choose Yourself:

% % Relaxation Terms:
T1 =   600;
T2 =    80;

% % RF pulses:
% rfa = [60,   90,  90]; % flip angles
% pha = [ 0,    0,   0]; % phase angles
% tau = [ 0,   16,  20]; % time before them
% rfa = [90, 180, 60, 90, 60]; % flip angles
% pha = [ 0,   0,  0,  0,  0]; % phase angles
% tau = [ 0,  10, 15, 15, 20]; % time before them

% % CPMG Weigel Paper: 
rfa = [90, 60, 60, 60, 60, 60]; % flip angles
pha = [90,  0,  0,  0,  0,  0]; % phase angles
tau = [ 0,  5, 10, 10, 10, 10]; % time before them
Npulses = length(rfa);       % how many
TR  = 1.5 * sum(tau); % total time of simulation

% % % % CPMG:
% % rfa = [90, 180, 180, 180, 180, 180]; % flip angles
% % pha = [ 0,   0,   0,   0,   0,   0]; % phase angles
% % tau = [ 0,   5,  10,  10,  10,  10]; % time before them
% % Npulses = length(rfa);       % how many
% % TR  = 1.5 * sum(tau); % total time of simulation

% % % % Carr-Purcell:
% % rfa = [90, 180, 180, 180, 180, 180]; % flip angles
% % pha = [ 0,  90,  90,  90,  90,  90]; % phase angles
% % tau = [ 0,   5,  10,  10,  10,  10]; % time before them
% % Npulses = length(rfa);       % how many
% % TR  = 1.5 * sum(tau); % total time of simulation

% % Spins: 
N   =  2000;
dw  =  12*pi/(2*pi*max(tau));
dws = linspace(-dw, dw, N); % create a population of spins with equally
                            % distributed off-resonances
%dw  = 0.01;                           
%dws = sqrt(dw).*randn(N,1); % create a population of spins with 
%                           % off-resonances drawn from a normal distribution

%% Create event matrix for plotting
dt  = 0.2; % how finely grained you want your simulation

% Events matrix with:
Tpoints = ceil(TR/dt);
% % % 1 = timepoint | 2 = RF.FA | 3 = RF.PhA
RFevents = zeros(3, Tpoints);
% Populate event matrix with the timpoints
RFevents(1, :) = linspace(0, TR, Tpoints);
% Populate eventmatrix with the rf pulses and phase angles
for i = 1:Npulses
	RFevents(2, floor(sum(tau(1:i))/dt)+1) =  rfa(i);
    RFevents(3, floor(sum(tau(1:i))/dt)+1) =  pha(i);
end

% Plot their evolution in time
M   = repmat([0 0 1]', 1, N);    % The m vectors
% separatingSpins = linspace(-0.5, 0.5, N);

% Figure
figure('Position', [10,10,1200,800])
nr = 5; nc = 11;
eventsPos = [ 1,  5]; % event matrix
signalPos = [ 7, 11]; % signal
simul3DPos = [  12, 49]; % the 3D view 
simul2DPos = [  18, 54]; % the 2D projection view

% Events
subplot(nr,nc,eventsPos)
scatter(RFevents(1,:), RFevents(3,:), 'ro', 'filled'), hold on % phase angles
stem(   RFevents(1,:), RFevents(2,:), 'b')
grid on
xlim([-1, TR])
ylabel('RF angle (deg)')
xlabel('t (ms)')


% Magnetic moments (the 3D plot)
subplot(nr,nc,simul3DPos)
c = linspace(1,2,N); % % Colors for the magnetic moments tips
axisl = 1.1;
createAxis(axisl);
axis square; %xlim([-1 1]); ylim([-1 1]); zlim([-1 1]); 
view(110, 10); %(110,10)
vecTip3D = scatter3(M(1, :), M(2, :), M(3, :), [], c, 'filled');

% Magnetic moments (the 2D plot)
subplot(nr,nc,simul2DPos)
plot([0 0], [1 -1], 'k--'), hold on, plot([1 -1], [0 0], 'k--')
vecTip2D = scatter(M(1, :), M(2, :), [], c, 'filled');
axis equal
xlim([-axisl axisl]); ylim([-axisl axisl])
xlabel('x') ; ylabel('y'); 
grid on

% Signal
subplot(nr,nc,signalPos)
% Add T2 exponential decay curve
plot(RFevents(1,:), exp(-RFevents(1,:)/T2)), hold on
% Add the RF pulses
stem(RFevents(1,:), RFevents(2,:)~=0, 'k--', 'marker', 'none')
grid on
xlim([-1, TR])
ylim([0, 1])
ylabel('|M_{xy}|')
xlabel('t (ms)')
legend('-e^{-t/T_2}')


pause


for i = 1:Tpoints
    
    % Delete history of vectors
    delete(vecTip3D); delete(vecTip2D);
    
    % Plot events
    subplot(nr,nc,eventsPos)
    event = stem(RFevents(1,i), 90, 'k');
    
    % DO RF PULSE
    if RFevents(2,i) ~= 0
        for j = 1:10                        
            Mtemp = Rotz( deg2rad(RFevents(3,i))) * ...
                     Rotx(-deg2rad(RFevents(2,i)*j/10)) * ...
                      Rotz(-deg2rad(RFevents(3,i))) * M;
                  
            % % % THE 3D plot
            subplot(nr,nc,simul3DPos)
            vecTip3D  = scatter3(Mtemp(1, :), Mtemp(2, :), Mtemp(3, :), ...
                           [], c, 'filled');
            hold on
                       
            % The 2D plot
            subplot(nr,nc,simul2DPos)
            vecTip2D = scatter(Mtemp(1, :), Mtemp(2, :), ...
                               [], c, 'filled');
            
            pause(0.1)
            delete(vecTip3D); delete(vecTip2D);
        end
        M = Rotz( deg2rad(RFevents(3,i))) * ...
                Rotx(-deg2rad(RFevents(2,i))) * ...
                    Rotz(-deg2rad(RFevents(3,i))) * M ;
    else
	% DEPHASE if not RF pulse events
        % Start dephasing for all magnetic moment vectors
        for j = 1:N
            M(:, j) = Drel(dt, T1, T2) * ...
                      Rotz(-2*pi*dt*dws(j)) * M(:, j) + ...
                      Drelz(dt, T1, 1);
        end
    end
    
    % % % The 3D plot
    subplot(nr,nc,simul3DPos)
    vecTip3D  = scatter3(M(1, :), M(2, :), M(3, :), ...
                         [], c, 'filled'); 
    hold on

    % The 2D plot
    subplot(nr,nc,simul2DPos)
    vecTip2D = scatter(M(1, :), M(2, :), [], c, 'filled');
    
    % Signal
    subplot(nr,nc,signalPos)
    scatter(RFevents(1,i), sqrt((sum(M(1,:))./N).^2 + ...
                                (sum(M(2,:))./N).^2), 'k.')

    hold on
    grid on
    xlim([0, TR])
    
    % delets
    pause(0.001)
    delete(event)
    
end

% % % end % FUNCTION END

