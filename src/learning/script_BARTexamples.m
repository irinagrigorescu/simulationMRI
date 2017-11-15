% % % % 
% % % % This script looks at BART reconstructions
% % % % 

% Prerequisites
addpath(genpath('~/OneDrive - University College London/Work/Simulation/simulationMRF/possumHelpers/functions/'))
addpath(genpath('~/Tools/bart-0.3.01/'))
setenv('TOOLBOX_PATH','/Users/irina/Tools/bart-0.3.01/')
addpath(genpath('~/Tools/MatlabStuff/'))

FOVx  = 256; FOVy = FOVx;
gamma = 2.68*1E08;
gammabar = gamma / (2*pi);

%% Look at radial reconsutrction with BART
figure
for i = 8:8
    
    % % Generate k-space trajectory with 
    % 4 * i radial spokes
    xSize = 500; ySize = 4 * i;
    scaleFactor = FOVx/xSize;
    traj_rad = bart(['traj -r -x', num2str(xSize), ' -y', num2str(ySize)]);

    % % Generate kspace data
    % increase the reconstructed FOV a bit
    traj_rad2 = bart(['scale ', num2str(scaleFactor)], traj_rad);

    % simulate 1-channel k-space data
    ksp_sim = bart('phantom -k -s1 -t', traj_rad2);

    % % Inverse gridding reconstruction
    % Apply inverse gridding
    igrid = bart('nufft -i -t', traj_rad2, ksp_sim);
    classicExample = abs(igrid);
    
    % Apply inverse gridding not taking into account the existence of
    % different spokes, but having all data together
    traj_rad3 = reshape(traj_rad2, [3, xSize*ySize]);
    ksp_sim3 = reshape(ksp_sim, [1, xSize*ySize]);
    igrid3 = bart('nufft -i -t', traj_rad3, ksp_sim3);
    myExample = abs(igrid3);
    
    % Start plotting
    subplot(2,2,1)
    scatter(traj_rad2(1, :), traj_rad2(2,:), 2, 'r'), 
    title(['Radial trajectory with ', num2str(ySize), ' spokes'])
    axis equal

    subplot(2,2,2)
    scatter(traj_rad2(1, :), traj_rad2(2,:), 1, log(ksp_sim(:))), 
    colormap(gray), title('Kspace data'), axis equal

    subplot(2,2,3) 
    imshow(myExample, [])
    title('Inverse Gridding'), axis equal

    subplot(2,2,4) 
    imshow(myExample-classicExample, [])
    title({'Difference between the', '2 reconstructions'})
    colorbar; axis equal 
    
    pause(0.1)
end

%% Look at spiral reconstruction with BART
% % Generate k-space trajectory with Nro points

alpha1 = 1 / (gammabar*0.3*230*1E-05); %1/(gammabar*L*n_theta*dt)
alpha2 = (2*pi) / (230*1E-05);         %(2*pi)/(n_theta*dt)

for i = 1:1
    % Number of readout points
    Nro = 6000.*i; 
    
    t = linspace(0, 6, Nro); 
    powerToRise = 1;
    traj_spiral = zeros(3, Nro); 
%     traj_spiral(1, :) = 0.53 .* (t.^powerToRise) .* cos(1340.*t);
%     traj_spiral(2, :) = 0.53 .* (t.^powerToRise) .* sin(1340.*t);
    traj_spiral(1, :) = alpha1 .* (t.^powerToRise) .* cos(alpha2.*t);
    traj_spiral(2, :) = alpha1 .* (t.^powerToRise) .* sin(alpha2.*t);

    % scaling factor to get to the desired FOV
    scaleFactor = FOVx/(2*max(max(max(abs(traj_spiral)))));

    % % Generate kspace data
    % increase the reconstructed FOV a bit
    traj_spiral2 = bart(['scale ', num2str(scaleFactor)], traj_spiral);

    % simulate 1-channel k-space data
    ksp_sim = bart('phantom -k -s1 -t', traj_spiral2);

    % % Inverse gridding reconstruction
    % Apply inverse gridding
    igrid = bart('nufft -i -t', traj_spiral2, ksp_sim);
    recoSpiral = abs(igrid);

    % % Start plotting
    figure
    subplot(2,2,1)
    scatter(traj_spiral2(1, :), traj_spiral2(2,:), 2, 'r'), 
    title(['Spiral trajectory with ', num2str(Nro), ' RO points'])
    axis equal

    subplot(2,2,2)
    scatter(traj_spiral2(1, :), traj_spiral2(2,:), 1, log(ksp_sim(:))), 
    colormap(gray), title('Kspace data'), axis equal

    subplot(2,2,3) 
    imshow(recoSpiral, [])
    title('Inverse Gridding'), axis equal
 
    pause(0.1)
end


%% Try to see if directly modifying the k-space with a mask affects it
% Shepp-logan phantom
P = phantom('Modified Shepp-Logan', FOVx/2);
% Shepp-logan fourier transformed
P_kspace = fftshift(fft2(P));

figure
for i = 10:-1:1
    % A circle mask to apply over it
    maskCircle = createCirclesMask(P, [FOVx/4 FOVx/4], 10.*i);
    % Kspace masked
    P_kspace_mask = P_kspace.*maskCircle;
    % Phantom inverse fourier transfromed
    P_ifft = ifft2(ifftshift(P_kspace_mask));

    % Plot phantom
    subplot(2,2,1)
    imagesc(P), colormap gray, axis off, axis equal
    % Plot k-space of phantom
    subplot(2,2,2)
    imagesc(abs(log(P_kspace))), colormap gray, axis off, axis equal
    % Plot reconstructed phantom
    subplot(2,2,3)
    imagesc(abs(P_ifft)), colormap gray, axis off, axis equal
    % Plot masked k-space
    subplot(2,2,4)
    imagesc(abs(log(P_kspace_mask))), colormap gray, axis off, axis equal
    
    pause
end