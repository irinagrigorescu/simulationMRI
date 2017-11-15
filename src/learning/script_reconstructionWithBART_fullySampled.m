% % % % 
% % % % This script handles a simulation made with JEMRIS
% % % % and tries to reconstruct with BART
% % % % 

% Set environment variable for BART
setenv('TOOLBOX_PATH','/Users/irina/Tools/bart-0.3.01/')

% Prerequisites
addpath(genpath('~/OneDrive - University College London/Work/Simulation/simulationMRF/possumHelpers/functions/'))
addpath(genpath('~/Tools/bart-0.3.01/'))
% Addpaths for EPG
addpath(genpath('~/OneDrive - University College London/Work/Simulation/simulationMRF/src/extendedPhaseGraph/'))
% Addpaths for possum helpers
addpath(genpath('~/OneDrive - University College London/Work/Simulation/simulationMRF/possumHelpers/functions/'))
% Addpaths for helpers
addpath(genpath('~/OneDrive - University College London/Work/Simulation/simulationMRF/src/helpers/'))
% Addpaths for algs
addpath(genpath('~/OneDrive - University College London/Work/Simulation/simulationMRF/src/algs/'))
% Addpaths for preprocess
addpath(genpath('~/OneDrive - University College London/Work/Simulation/simulationMRF/src/preprocess/'))
% Addpaths for perlin
addpath(genpath('~/OneDrive - University College London/Work/Simulation/simulationMRF/src/perlin/'))

%%

% JEMRIS simulation
folderName = '~/jemrisSims/spiral/testFullySampled/';%testFullySampled/';

% Set number of repetitions for the fully sampled
NReps = 1;
% Set number of redout points per TR
nROPoints = 6000; %125; %6000;
% Set number of TR blocks
NTRs  = 5; %20; %5;
% Set effective image size
Nx   = 64; Ny = Nx; 

% Allocate memory for k-space coordinates in x,y,z directions
kspace_x = zeros(nROPoints, NReps, NTRs); % NROpoints x Nrepetitions x NTRs
kspace_y = zeros(nROPoints, NReps, NTRs);
kspace_z = zeros(nROPoints*NReps, 1);

for i = 1:NReps
    % Get K-space coordinates from JEMRIS simulation
    % along with other sequence information
    filenameSequence = [folderName, 'seq', num2str(i),'.h5']; 
    t   = h5read(filenameSequence,'/seqdiag/T');    % Timepoints
    RXP = h5read(filenameSequence,'/seqdiag/RXP');  % Read-Out (-1 no RO / 0 RO)
    TXM = h5read(filenameSequence,'/seqdiag/TXM');  % RF pulse FA  in rad
    TXP = h5read(filenameSequence,'/seqdiag/TXP');  % RF pulse PhA in rad
    GX  = h5read(filenameSequence,'/seqdiag/GX');   % GX (mT/m)
    GY  = h5read(filenameSequence,'/seqdiag/GY');   % GY (mT/m)
    GZ  = h5read(filenameSequence,'/seqdiag/GZ');   % GZ (mT/m)
    KX  = h5read(filenameSequence,'/seqdiag/KX');   % KX (rad/mm)
    KY  = h5read(filenameSequence,'/seqdiag/KY');   % KY (rad/mm)
    KZ  = h5read(filenameSequence,'/seqdiag/KZ');   % KZ (rad/mm)
    B   = h5read(filenameSequence,'/seqdiag/META');
    pulseSequenceJ = [t./1000, TXM      , zeros(size(TXM)), ...
                                          zeros(size(TXM)), ...
                      RXP+1   , GX./1000, ...
                                GY./1000, ...
                                GZ./1000];

    % Calculate number of readout points - this will give you 300
    if nROPoints ~= (size(RXP(RXP==0),1)./NTRs)
        disp('Warning - unexpected number of READOUT points');
    end

    % Retrieve all k-space coordinates where READOUT has happened
    % for each sequence, meaning it will be 300 x 5 (NROPoints x NTRs)
    kspx = KX(RXP==0); kspy = KY(RXP==0);

    % For each TR block calculate the kspace coordinates by bringing them back
    % to 1/cm from rad/mm
    % My signals will be in the form of NRO x NReps x NTRs 
    % for example, my signals will be   300 x  48   x  5  
    for j = 1:NTRs
        % NROpoints x Nrepetitions x NTRs
        kspace_x(:,i,j) = (kspx((j-1)*nROPoints+1 : j*nROPoints) ...
                           ./ (2*pi)) .* 1E+01; %cm^-1
        kspace_y(:,i,j) = (kspy((j-1)*nROPoints+1 : j*nROPoints) ...
                           ./ (2*pi)) .* 1E+01; %cm^-1
    end
    
end

% Reshape to fully sampled data
kspace_x = reshape(kspace_x, [nROPoints*NReps, NTRs]);
kspace_y = reshape(kspace_y, [nROPoints*NReps, NTRs]);

% Plot them for each spiral turn
figure(111)
for i = 1:NTRs
    % Plot the coordinates before scaling
    subplot(1,2,1); xl = 4.5;
    scatter(kspace_x(:, i),kspace_y(:, i),'.'), axis equal, axis square;
    xlabel('k_x cm^{-1}'); ylabel('k_y cm^{-1}');
    title(['Spiral ', num2str(i)])
    xlim([-xl xl]); ylim([-xl xl]);

    % Transform the kspace coordinates to the desired effective imaging matrix
    % by scaling the values so that they fit into the matrix
    dist_ksp = sqrt(   (kspace_x(end,i) - kspace_x(1,i)).^2 + ...
                       (kspace_y(end,i) - kspace_y(1,i)).^2  );
    kspace_x(:,i) = kspace_x(:,i) .* Nx./(2*dist_ksp);
    kspace_y(:,i) = kspace_y(:,i) .* Ny./(2*dist_ksp);

    % Plot the coordinates after scaling
    subplot(1,2,2)
    scatter(kspace_x(:, i), kspace_y(:, i),'.'), axis equal, axis square;
    xlabel('N_x'); ylabel('N_y');
    title(['Spiral ', num2str(i)])
    xlim([-Nx/2 Nx/2]); ylim([-Ny/2 Ny/2]);

    drawnow
end


%%  GET THE SIGNAL VALUES
% Allocate memory for kspace values
kspaces = zeros(nROPoints, NReps, NTRs); % NROpoints x Nrepetitions x NTRs

for i = 1:NReps
    % Get signal values from JEMRIS simulation
    filenameSignal = [folderName, 'signals', num2str(i),'.h5']; 

    % % % Read in the vector values
    A = (h5read(filenameSignal, '/signal/channels/00'))';
    % % % Get timings of signal readout
    timingRO = h5read(filenameSignal, '/signal/times');
    % % % Readout points
    Nro = size(timingRO, 1);
    % % % Sort them and take signal in sorted order
    [timingRO,J] = sort(timingRO); Mvecs = A(:,:);
    % % % Calculate timings between each line of k-space and
    % % % between different slices
    d = diff(diff(timingRO));
    d(d<1e-5) = 0;

    % % % And store them in the solumn vector of indices I
    % % % together with index 0 and last index
    I = [0; find(d) + 1; length(timingRO)];

    % Get kspace as one signal
    for j = 1:NTRs
        kspaces(:,i,j) = Mvecs((j-1)*nROPoints+1 : j*nROPoints, 1).' + ...
              sqrt(-1) * Mvecs((j-1)*nROPoints+1 : j*nROPoints, 2).';

    end
end

% Reshape kspace data to have all the readout points from all repetitions
kspaces = reshape(kspaces, [nROPoints*NReps, NTRs]);


%% Let's try BART
% inverse gridding

% Save signal from image space from one voxel from the image
signalImage = zeros(NTRs, 1);

figure(223)
for i = 1:NTRs
    tic
    % Retrieve current values for current reconstruction
    currentKspace_x = squeeze(kspace_x(:, i));
    currentKspace_y = squeeze(kspace_y(:, i));
    currentKspaceToReco = squeeze(kspaces(:, i)).';
    
    % Prepare trajectory for current image
    traj_kspace = [ currentKspace_x, ...
                    currentKspace_y, ...
                    kspace_z               ].';
                
    % Reconstruct with BART
    igrid = bart('nufft -i -t', traj_kspace, currentKspaceToReco);
    toc
    
    % Save value of signal at a randomly chosen point in k-space
    signalImage(i) = igrid(Nx/2, Ny/2); %igrid(32, 32); %Nx/2, Ny/2));
    
    imagesc(abs(igrid))
    title(['Spiral reconstruction ', num2str(i)])
    axis square; axis off;
    colormap gray;
    colorbar;
    
    pause%(0.1)
end


%%
% % % % % % % % % % % Load EPG Dictionary
% Load data from dictionary file 
varToLoad = load(['/Users/irina/OneDrive - University College London/', ...
                  'Work/Simulation/simulationMRF/src/', ...'
                  'data_EPGDictionary/', ...
                  'T11000T2100NTR500Mariya.mat']);
% Get values from loaded variable
F0states = varToLoad.F0states(2:end); % 2:end because first is after IR 
% Material tuples
materialTuplesEPG = varToLoad.materialTuples;

% Create 2 channel data (real + imag)
F0statesReIm = [ real(F0states), imag(F0states)];

% Normalize signals
[dictionaryMRFFISPNorm2Ch, dictionaryMRFFISPNorm] = ...
        normalizeSignals(F0statesReIm);

% Calculate absolute signal
sigMRFFISPAbsNorm = abs(dictionaryMRFFISPNorm(:, 1:end));



%% Manipulate the saved signal
[~, signalImage2] = normalizeSignals([real(signalImage.'), imag(signalImage.')]);

figure
subplot(2,1,1)
plot(1:NTRs, abs(signalImage2), 'o-')
hold on
plot(1:NTRs, abs(dictionaryMRFFISPNorm(1:NTRs)), 'o-')
xlabel('TR block number');
ylabel('normalised signals (a.u.)');
legend('jemris', 'epg')

subplot(2,1,2)
plot(1:NTRs, abs(signalImage2) - abs(dictionaryMRFFISPNorm(1:NTRs)), 'o-')
xlabel('TR block number');
ylabel('Jemris - Dictionary difference');








