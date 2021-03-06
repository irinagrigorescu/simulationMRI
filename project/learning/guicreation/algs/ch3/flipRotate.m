function [ mu, muz, muxy, b1field ] = flipRotate( gamma, B0, B1, ...
                    phaseAngle, dw, tau, ...
                    rotatingFrame, rotateB1 )
% % % % IRINA GRIGORESCU
% % % % DATE: 28-Nov-2016
% % % % CHAPTER 3.3.1
% % % % 
% % % % This function calculates the trajectory of a precessing spin
% % % % INPUT PARAMETERS
% % % %     gamma = 
% % % % 
% % % %     B0 = 
% % % % 
% % % %     B1 = 
% % % % 
% % % %     phaseAngle = will be the phase angle
% % % % 
% % % %     dw = off-resonance frequency
% % % % 
% % % %     tau = 
% % % % 
% % % %     rotatingFrame = 
% % % % 
% % % %     rotateB1 = 
% % % % 
% % % % 
% % % % OUTPUT PARAMETERS: - Trajectory of spin through time
% % % %     mu = magnetic moment vector trajectory
% % % % 
% % % %     muz = mu projection on z axis trajectory 
% % % % 
% % % %     muxy = mu projection on transverse plane trajectory
% % % % 
% % % %     b1field = b1 field trajectory
% % % % 

%% Rotation Matrices
% % % % Creating some local functions for rotations of vectors

% % % % Rotation about x
Rx = @(omega,t) [ 1         0            0       ; 
                  0   cos(omega*t)  -sin(omega*t) ;
                  0   sin(omega*t)   cos(omega*t) ];
% % % % Rotation about y
Ry = @(omega,t) [  cos(omega*t)  0   sin(omega*t) ; 
                        0        1         0      ;
                  -sin(omega*t)  0   cos(omega*t) ];
% % % % Rotation about z
Rz = @(omega,t) [ cos(omega*t) -sin(omega*t) 0 ; 
                  sin(omega*t)  cos(omega*t) 0 ;
                         0            0      1 ];              

%% Useful parameters

% Unit vectors:
x = [1 0 0]'; y = [0 1 0]'; z = [0 0 1]';
% Timestep for update
dt = 1e-6; 
% Number of time steps 
N = ceil(tau/dt);
% Precession angular frequency
omega0 = gamma * B0;
% Spin precession frequency generated by the circularly polarized rf field
omega1 = gamma * B1;

                    
%% Calculated Variables
% % % % As output we will have vectors of 3xN dimension corresponding to
% % % % the 3D coordinates for each timestep from 1 to N

% Initialise variables
mu      = zeros(3,N);  % magnetic moment
muz     = zeros(3,N);  % projection of magnetic moment on z
muxy    = zeros(3,N);  % projection of magnetic moment on xy
b1field = zeros(3,N);  % b1 field is a rotating one

mu(:,1)       = 0*x + 0*y + 1*z;
muz(:,1)      = mu(:,1);         %0*x + 0*y + 1*z;
muxy(:,1)     = 0*x + 0*y + 0*z;
b1field(:,1)  = 1*x + 0*y + 0*z;

%% Calculate the incremental phase angle of the rf pulse
ph = phaseAngle ./ tau;

%% Create rotation functions 

% Rz(ph, tau) makes an instantaneous rotation about z with the phase angle
% ph; rotations are clockwise
switch rotateB1
    case '+x'
        b1field(:,1)  = Rz(ph, tau) * (1*x + 0*y + 0*z);
        RotB1 = @(omega, t) Rz(ph, tau) * Rx(omega, t) * Rz(-ph, tau);
    case '-x'
        b1field(:,1)  = Rz(ph, tau) * (-1*x + 0*y + 0*z);
        RotB1 = @(omega, t) Rz(ph, tau) * Rx(-omega, t) * Rz(-ph, tau);
    case '+y'
        b1field(:,1)  = Rz(ph, tau) * (0*x + 1*y + 0*z);
        RotB1 = @(omega, t) Rz(ph, tau) * Ry(omega, t) * Rz(-ph, tau);
    case '-y'
        b1field(:,1)  = Rz(ph, tau) * (0*x - 1*y + 0*z);
        RotB1 = @(omega, t) Rz(ph, tau) * Ry(-omega, t) * Rz(-ph, tau);
end

switch rotatingFrame
    % When not in the rotating frame things go crazy
    % The rotationg about z with dw is due to local inhomogeneities
    case 'no'
        RotAxis = @(omega, ph, dw, t) Rz(dw, t) * Rz(-omega,t);

    % When you are in the rotating frame things look stationaryti
    % and rotation is caused only by off-resonance stuff
    case 'yes'
        RotAxis = @(omega, ph, dw, t) Rz(dw, t);
end


%% Calculate vector movement history
N
for i = 2:N
    % Update magnetic moment vector position for each timestep
    mu(:,i)   = RotAxis(-omega0, ph, dw, i*dt) * ...
                RotB1(-omega1, i*dt) * mu(:,1);
    % Update magnetic moment vector projections for each timestep
    muz(3,i)  = mu(3,i); 
    muxy(1:2,i) = mu(1:2,i);
    % Update B1 field vector position for each timestep
    b1field(:,i) = RotAxis(-omega0, 0, 0, i*dt) * b1field(:,1);
end


end

