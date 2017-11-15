% % % % 
% % % % This script handles a simulation made with JEMRIS
% % % % and tries to reconstruct with BART
% % % % 

% Prerequisites
addpath(genpath('~/OneDrive - University College London/Work/Simulation/simulationMRF/possumHelpers/functions/'))
addpath(genpath('~/Tools/bart-0.3.01/'))
setenv('TOOLBOX_PATH','/Users/irina/Tools/bart-0.3.01/')

% JEMRIS simulation
folderName = '~/jemrisSims/spiral/';

% Get K-space coordinates from JEMRIS simulation
% along with other sequence information
filenameSequence = [folderName, 'seq.h5']; 
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

% Set effective image size
Nx   = 256; Ny = Nx; 
% Set number of TR blocks
NTRs = 1;
% Calculate number of readout points (equal number of RO points per TR
% block)
nROPoints = size(RXP(RXP==0),1)./NTRs;

% Retrieve all k-space coordinates where READOUT has happened
kspx = KX(RXP==0); kspy = KY(RXP==0);

% Allocate memory for k-space coordinates in x,y,z directions
kspace_x = zeros(nROPoints, NTRs);
kspace_y = zeros(nROPoints, NTRs);
kspace_z = zeros(nROPoints, 1);

% For each TR block calculate the kspace coordinates by bringing them back
% to 1/cm from rad/mm
for i = 1:NTRs
    kspace_x(:,i) = (kspx((i-1)*nROPoints+1:i*nROPoints) ./ (2*pi)) ...
                                                         .* 1E+01; %cm^-1
    kspace_y(:,i) = (kspy((i-1)*nROPoints+1:i*nROPoints) ./ (2*pi)) ...
                                                         .* 1E+01; %cm^-1
end

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


%%
% Get signal values from JEMRIS simulation
filenameSignal = [folderName, 'signals.h5']; 

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
kspaces = zeros(NTRs, nROPoints);
for i = 1:NTRs
    kspaces(i,:) = Mvecs((i-1)*nROPoints+1:i*nROPoints, 1).' + ...
        sqrt(-1) * Mvecs((i-1)*nROPoints+1:i*nROPoints, 2).';
end


% Get kspace in matrix directly
kspacesAsImage = [];
for i = 1:length(I)-1
    J = [I(i)+1:I(i+1)]';

    kspacesAsImage(:,i) = Mvecs(J, 1) + sqrt(-1) * Mvecs(J, 2);
end

figure, hold on
scatter(timingRO, Mvecs(:,1), '.')
scatter(timingRO, Mvecs(:,2), '.')
scatter(timingRO, Mvecs(:,3), '.')
xlabel('ms'); ylabel('M_i');
legend('M_x','M_y','M_z')


%%
% Let's try BART
% inverse gridding
figure(223)
for i = 1:1
    tic
    % Prepare trajectory for current image
    traj_kspace = [ reshape(kspace_x, [nROPoints*NTRs, 1]), ...
                    reshape(kspace_y, [nROPoints*NTRs, 1]), ... 
                    zeros(nROPoints*NTRs, 1)              ].';
                
    kspacesAll = reshape(kspaces, [nROPoints*NTRs, 1]).';
    % Reconstruct with BART
    igrid = bart('nufft -i -t', traj_kspace, kspacesAll);
    toc
    
    imagesc(abs(igrid))%./max(max(igrid))))
    title(['Spiral reconstruction ', num2str(i)])
    axis square; axis off;
    colormap gray;
    colorbar;
    
    pause(1)
end


 










