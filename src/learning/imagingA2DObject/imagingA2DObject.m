% function imagingA2DObject()
% 
% This is a script to simulate the GRE acquisition of an object
% The object is represented by position in x and y, T1 and T2
% 

addpath(genpath('~/Tools/MatlabStuff/'))
addpath(genpath('~/OneDrive - University College London/Work/Simulation/simulationMRF/src/helpers/'))

% % % % % % % % % % % % % % % 
% % % Prerequisites:
gamma    = 2.68*1E08;      %(rad/s/T)
gammabar = gamma./(2*pi);  %(1/s/T)
% Choose relaxation times:
T1 = [950, 600, 4500, 10000000000]./1E+03; % GM,WM,CSF
T2 = [100,  80, 2200, 10000000000]./1E+03;

% % % % % % % % % % % % % % % 
% % % Object
FOVx = 64; FOVy = 64;
[object2D, masks] = createObject(FOVx, FOVy, T1, T2, 0);

% % % % % % % % % % % % % % % 
% % % Imaging parameters:
dx         = 1*1E-03;       dy         = 1*1E-03;       % image voxel dim
Nx         = 64;            Ny         = 64;            % resolution of img
kspace     = zeros(Nx, Ny);                             % allocating space

[dkx, dky] = createImagingParameters(Nx, Ny, dx, dy);   % calculate 
                                                        % sampling
                                                        % distances in
                                                        % k-space
% % % % % % % % % % % % % % % 
% % % Sequence:
samplingFrequency = 1000000;
[dtSampling, Tsampling, dtG, dtPhase, ...
          GxAmplitude, GxArea, ...
          GxPreAmplitude, GxPreArea, ...
          GyAmplitude, dGyAmplitude] = createSequenceParametersGRE( ...
                    samplingFrequency, Nx, Ny, dkx, dky, gammabar);

%%
% % % Sequence events:
FA = 90; TE = 0.10; TR = 50;
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
        
        % Re-initialise the event matrix for each voxel
        eventMatrixPerVoxel = eventMatrix;
        
        % For each phase encoding line
        for phLine = Ny:-1:1 %Ny:-1:1
            % Initialise magnetisation vector 
            M = zeros(3, 1); M(3,1) = 1;
            
            % 1. Do RF pulse
            M = Rotx(-eventMatrixPerVoxel(1,2)) * M; 
            
            % 2. Relax until first gradient
            M = Drel(eventMatrixPerVoxel(2,1), ...
                chosenVoxel(3), chosenVoxel(4)) * M + ...
                Drelz(eventMatrixPerVoxel(2,1), chosenVoxel(3), 1);
            
            % 3. Prewinder x gradient and y gradient effects
            % on x
            phiGradX = gamma * eventMatrixPerVoxel(3,4) * ...
                       chosenVoxel(1) * ...
                       eventMatrixPerVoxel(3,1); % gamma*Gx*posx*dt
            % on y
            phiGradY = gamma * eventMatrixPerVoxel(3,5) * ...
                       chosenVoxel(2) * ...
                       eventMatrixPerVoxel(3,1); % gamma*Gy*posy*dt
            % do rotation
            M = Rotz(-phiGradY) * Rotz(-phiGradX) * M;
            % Relax for this time
            M = Drel(eventMatrixPerVoxel(3,1), ...
                     chosenVoxel(3), chosenVoxel(4))  * M + ...
                Drelz(eventMatrixPerVoxel(3,1), chosenVoxel(3), 1);
            
            % 4. Readout
            for roPoint = 1:Nx
                % gamma*Gx*posx*dt
                phiGradRO = gamma * eventMatrixPerVoxel(3+roPoint,4) * ...
                         chosenVoxel(1) * eventMatrixPerVoxel(3+roPoint,1); 
                % Apply gradient effect
                M = Rotz(-phiGradRO) * M;
                % Relax for this time
                M = Drel(eventMatrixPerVoxel(3+roPoint,1), ...
                         chosenVoxel(3),   ...
                         chosenVoxel(4)   ) * M + ...
                    Drelz(eventMatrixPerVoxel(3+roPoint,1), ...
                          chosenVoxel(3), 1);
                % Do Readout
                if (eventMatrixPerVoxel(3+roPoint,3) == 1)
                    mCplx = M(1) + 1i.*M(2);
                    kspace(phLine, roPoint) = ...
                                    kspace(phLine, roPoint) + mCplx;
                end
                                
            end
            
            % 5. Relaxation after Readout - not necessary

            % 6. Increase Gradient y
            eventMatrixPerVoxel(3,5) = eventMatrixPerVoxel(3,5) + ...
                                       dGyAmplitude;
            
        end
        
    end 
    
    t = toc;
    fprintf('Done with line %d in t = %d seconds\n', voxX, t);
end

%%
kspaceUnmodified = kspace;
figure, imagesc(abs(kspace))


%%
myImage = fftshift(ifft2(ifftshift(kspaceUnmodified)));

figure
subplot(1,2,1)
imagesc(abs((myImage))); 
axis square; 
colorbar
title('my kspace')

subplot(1,2,2)
imagesc(abs(((fftshift(fft2(squeeze(masks.mask2.*T1(2))))))))
axis square
colorbar

