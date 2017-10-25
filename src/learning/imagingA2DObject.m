% function imagingA2DObject()
% 
% This is a script to simulate the GRE acquisition of an object
% 
% 

addpath(genpath('~/Tools/MatlabStuff/'))
addpath(genpath('~/OneDrive - University College London/Work/Simulation/simulationMRF/src/helpers/'))

% % % % % % % % % % % % % % % 
% % % Prerequisites:
gamma    = 2.68*1E08;      %(rad/s/T)
gammabar = gamma./(2*pi);  %(1/s/T)

% % % % % % % % % % % % % % % 
% % % Object
FOVx = 64;  % size of the object in x
FOVy = 64;  % size of the object in y
posx = (-FOVx/2)/1E03 : 1/1E03 : (FOVx/2-1)/1E03; % x positions in meters
posy = (-FOVy/2)/1E03 : 1/1E03 : (FOVy/2-1)/1E03; % y positions in meters
[Xcoord, Ycoord] = meshgrid(posx, posy); Ycoord = Ycoord.'; % coordinates

% My object has 4 parameters: xcoord, ycoord, T1, T2
object2D = zeros(FOVx, FOVy, 4); 
object2D(:,:,1) = Xcoord;
object2D(:,:,2) = Ycoord;

% Create masks for 3 objects in the image
mask1 = createCirclesMask([FOVx FOVy], ...
                          [FOVx/2-FOVx/4, FOVy/2-FOVy/4], 8);
                                                 % upper left circle
mask2 = createCirclesMask([FOVx FOVy], ...
                          [FOVx/2, FOVy/2], 8); % centre circle
mask3 = createCirclesMask([FOVx FOVy], ...
                          [FOVx/2+FOVx/4, FOVy/2+FOVy/4], 8);
                                                 % lower right circle

% Choose relaxation times:
T1 = [950, 10000000000, 4500]./1E+03; % GM,WM,CSF
T2 = [100, 10000000000, 2200]./1E+03;

% Create object
object2D(:,:,3) = mask1.*T1(2);
object2D(:,:,4) = mask1.*T2(2);

% Figure:
figure
subplot(1,2,1)
imagesc(squeeze(object2D(:,:,3)));
axis off; axis square; colorbar
title('T_1')

subplot(1,2,2)
imagesc(squeeze(object2D(:,:,4)))
axis off; axis square; colorbar
title('T_2')

% % % % % % % % % % % % % % % 
% % % Imaging parameters:
Nx  = 64;        Ny  = 64;
dx  = 1*1E-03;   dy  = 1*1E-03;
dkx = 1/(dx*Nx); dky = 1/(dy*Ny);
kspaceReal = zeros(Nx, Ny);
kspaceImag = zeros(Nx, Ny);

% % % % % % % % % % % % % % % 
% % % Sequence:
Tsampling  = 0.004;          % sampling time per line in seconds
dtSampling = Tsampling / Nx; % dt sampling
dtG        = Tsampling/2;    % duration of each square gradient
dtPhase    = dtG;     % dt phase encoding

GxAmplitude    = dkx/(gammabar*dtSampling);   % T/m
GxArea         = GxAmplitude * Tsampling;     % Ts/m
GxPreArea      = GxArea / 2;                  % Ts/m
GxPreAmplitude = GxPreArea / dtG;             % T/m
dGyAmplitude   = dky/(gammabar*dtPhase);      % T/m
GyAmplitude    = (Ny/2-1)*dGyAmplitude;       % T/m

% % % Sequence events:
FA = 90; TE = 0.012; TR = 50;
eventMatrix = zeros(3+Nx, 5); % 1=dTime, 2=RFpulse, 3=RO, 4=Gx, 5=Gy
% 1. RF pulse
eventMatrix(1,2) = deg2rad(FA);
% 2. Waiting Time
eventMatrix(2,1) = TE - Tsampling; % dt (seconds)
% 3. Gx prewinder and Gy
eventMatrix(3,1) =  dtG;            % dtG
eventMatrix(3,4) = -GxPreAmplitude; %  Gx prewinder
eventMatrix(3,5) = -GyAmplitude;    %  Gy 
% 4. Gx readout
for roPoint = 1:Nx
    eventMatrix(3+roPoint,1) =   dtSampling;   % Tsampling/Nx
    eventMatrix(3+roPoint,3) =            1;   % readout
    eventMatrix(3+roPoint,4) =  GxAmplitude;   % Gx
end
% 5. Relaxation after Readout
eventMatrix(3+roPoint+1,1)   =  TR - TE - Tsampling/2;   % seconds

% figure; hold on
% plot(eventMatrix(1:20,1), eventMatrix(1:20,4), '*'); 
% plot(eventMatrix(1:20,1), eventMatrix(1:20,5), 'o'); 

% % % % % % % % % % % % % % % 
% % % Simulate:

% for each voxel in the 2D object
for voxX = 1:FOVx %FOVx/4 % 1:FOVx
    tic;
    for voxY = 1:FOVy % 1:FOVy
        
        % Choose voxel from the initial image
        chosenVoxel = squeeze(object2D(voxX, voxY, :));
        % If it's zero move to the next one
        if chosenVoxel(3) == 0 && chosenVoxel(4) == 0
            continue
        end
        
        % For each phase encoding line
        for phLine = Ny:-1:1
            % Initialise magnetisation vector 
            M = zeros(3, 1); M(3,1) = 1;
            
            % 1. Do RF pulse
            M = Rotx(-eventMatrix(1,2)) * M; % rf pulse on first position
            % 2. Relax until first gradient
            M = Drel(eventMatrix(2,1), chosenVoxel(3), chosenVoxel(4)) ...
                     * M + ...
                Drelz(eventMatrix(2,1), chosenVoxel(3), 1);
            % 3. Prewinder x gradient and y gradient effects
            % on x
            phiGradX = gamma * eventMatrix(3,4) * ...
                       chosenVoxel(1) * eventMatrix(3,1); % gamma*Gx*posx*dt
            % on y
            phiGradY = gamma * eventMatrix(3,5) * ...
                       chosenVoxel(2) * eventMatrix(3,1); % gamma*Gy*posy*dt
            M = Rotz(-phiGradY) * Rotz(-phiGradX) * M;
            % Relax for this time
            M = Drel(eventMatrix(3,1), chosenVoxel(3), chosenVoxel(4)) ...
                     * M + ...
                Drelz(eventMatrix(3,1), chosenVoxel(3), 1);
            
            % 4. Readout
            for roPoint = 1:Nx
                % gamma*Gx*posx*dt
                phiGradRO = gamma * eventMatrix(3+roPoint,4) * ...
                            chosenVoxel(1) * eventMatrix(3+roPoint,1); 
                % Apply gradient effect
                M = Rotz(-phiGradRO) * M;
                % Relax for this time
                M = Drel(eventMatrix(3+roPoint,1), ...
                         chosenVoxel(3),   ...
                         chosenVoxel(4)   ) * M + ...
                    Drelz(eventMatrix(3+roPoint,1), chosenVoxel(3), 1);
                % Do Readout
                if (eventMatrix(3+roPoint,3) == 1)
                    m00 = abs(M(1) + 1i.*M(2));
                    kspaceReal(phLine, roPoint) = ...
                                    kspaceReal(phLine, roPoint) + ...
                                    m00 * cos(atan2(M(2),M(1)));
                    kspaceImag(phLine, roPoint) = ...
                                    kspaceImag(phLine, roPoint) + ...
                                    m00 * sin(atan2(M(2),M(1)));

                end
            end
            
            % 5. Relaxation after Readout - not necessary

            % 6. Increase Gradient y
            eventMatrix(3,5) = eventMatrix(3,5) + dGyAmplitude;
        end
        
    end
    
    t = toc;
    fprintf('Done with line %d in t = %d seconds\n', voxX, t);
end
%%
kspaceUnmodified = kspaceReal + 1i.*kspaceImag;
figure, imagesc(abs(kspaceImag))


%%
myImage = fftshift(ifft2(ifftshift(kspaceUnmodified)));

figure
subplot(1,2,1)
imagesc(abs(kspaceUnmodified)); 
axis square; 
colorbar
title('my kspace')

subplot(1,2,2)
imagesc(abs(((fftshift(fft2(squeeze(mask2.*T1(2))))))))
axis square
colorbar