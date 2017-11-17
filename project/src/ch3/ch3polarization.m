% % % % IRINA GRIGORESCU
% % % % DATE: 24-Nov-2016
% % % % CHAPTER 3.2.1
% % % % Polarization

addpath ../helpers

%% Parameters 
gamma = 2.68*(10^8);  % gyromagnetic ratio
B0 = 3;     %1e-10; % external magnetic fields in Tesla
omega0 = gamma .* B0; % rad/s
omega1 = omega0*1e-01;

%% Rotations matrices
Rz = @(omega,t) [ cos(omega*t) sin(omega*t) 0 ; 
                 -sin(omega*t) cos(omega*t) 0 ;
                        0            0      1 ];
                   
Rx = @(omega,t) [ 1         0            0       ; 
                  0   cos(omega*t)  sin(omega*t) ;
                  0  -sin(omega*t)  cos(omega*t) ];

Ry = @(omega,t) [ cos(omega*t)  0   sin(omega*t) ; 
                        0       1         0      ;
                 -sin(omega*t)  0   cos(omega*t) ];

%% Parameters to be calculated
angleOfRot = 90;
angle = angleOfRot * pi/180; % rad = degrees*pi/180
tau = angle / omega1; % deltaTheta = gamma * B1 * tau = omega1 * tau
dt = 1e-10;        % timestep for update
N = ceil(tau/dt); % number of time steps 
mu      = zeros(3,N);  % vector
muz     = zeros(3,N);  % vector
muxy    = zeros(3,N);  % vector
b1field = zeros(3,N);  % vector
xprime  = zeros(3,N);
yprime  = zeros(3,N);

% Unit vectors:
x = [1 0 0]'; y = [0 1 0]'; z = [0 0 1]';
% Magnetic moment:
mu(:,1)       = 0*x + 0*y + 1*z;
% Magnetic moment projections
muz(:,1)      = 0*x + 0*y + 1*z;
muxy(:,1)     = 0*x + 0*y + 0*z;
% B1 vector
b1field(:,1)  = 0*x + 0*y + 0*z;

%% Create rotation functions 
rotatingFrame = 'no';
rotateB1 = 'x';

switch rotateB1
    case 'x'
        b1field(:,1)  = 1*x + 0*y + 0*z;
        RotB1 = @(omega, t) Rx(omega, t);
        getMuxy = @(vec) [vec(1,2) vec(3,2) vec(2,2)]';
    case 'y'
        b1field(:,1)  = 0*x + 1*y + 0*z;
        RotB1 = @(omega, t) Ry(omega, t);
        getMuxy = @(vec) [vec(3,2) vec(1,2) vec(2,2)]';
end

switch rotatingFrame
    % When not in the rotating frame things go crazy
    case 'no'
        RotAxis = @(omega, t) Rz(-omega,t);

    % When you are in the rotating frame things look stationary
    case 'yes'
        RotAxis = @(omega, t) Rz(omega,t) * Rz(-omega,t);
end


%% Calculate vector movement history
for i = 2:N
    % Update magnetic moment vector position for each timestep
    mu(:,i)   = RotAxis(omega0, i*dt) * RotB1(omega1, i*dt) * mu(:,1);
    % Update magnetic moment vector projections for each timestep
    muz(:,i)  = cos(omega1*i*dt) .* muz(:,1);
    muxy(:,i) = sin(omega1*i*dt) .* (RotAxis(omega0, i*dt) * getMuxy(mu));
    % Update B1 field vector position for each timestep
    b1field(:,i) = RotAxis(omega0, i*dt) * b1field(:,1); 
end

%% Plot rotating frame
% Create static coordinate system
[xax, yax, zax] = createAxis();
% Viewer position
view([45 45 10]);

%%
% Draw for each timestep
for i = 2:N
    % Vector from (0,0,0) to (x,y,z) of updated position
    muvec = [[0 0 0]'  mu(:,i)]; 
    muxyvec = [[0 0 0]'  muxy(:,i)]; 
    muzvec = [[0 0 0]'  muz(:,i)]; 
    b1vec = [[0 0 0]'  b1field(:,i)]; 
    xprimevec = [[0 0 0]'  xprime(:,i)]; 
    yprimevec = [[0 0 0]'  yprime(:,i)]; 
    
    % Plot vector and arrow at the end
    % for mu vector
    muvecplot = plot3(muvec(1,:), muvec(2,:), muvec(3,:), 'r'); hold on
    mutipplot = scatter3(mu(1,i), mu(2,i), mu(3,i), 'r^'); hold on
    % for xy projection of mu vector
    muxyplot = plot3(muxyvec(1,:), muxyvec(2,:), muxyvec(3,:), '--r'); hold on
    muxytipplot = scatter3(muxy(1,i), muxy(2,i), muxy(3,i), 'r*'); hold on
    % for z projection of mu vector
    muzplot = plot3(muzvec(1,:), muzvec(2,:), muzvec(3,:), '--r'); hold on
    muztipplot = scatter3(muz(1,i), muz(2,i), muz(3,i), 'r*'); hold on
    % for B1 field
    b1vecplot = plot3(b1vec(1,:), b1vec(2,:), b1vec(3,:), '--k'); hold on
    b1tipplot = scatter3(b1field(1,i), b1field(2,i), b1field(3,i), 'ko'); hold on
    
    % Draw and pause
    drawnow
    pause(0.001)
    
    % Delete vector and arrow so that each update looks noew
    if i ~= N
        deletePlots(muvecplot, mutipplot, ...
                    muxyplot, muxytipplot, ...
                    muzplot, muztipplot, ...
                    b1vecplot, b1tipplot); 
    end
    
    % Plot history
    scatter3(b1field(1,i), b1field(2,i), b1field(3,i), 'k.'); hold on
end





