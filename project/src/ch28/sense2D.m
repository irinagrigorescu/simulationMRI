function [ reconstruction ] = sense2D( aliasedImages, sensProfiles, Rx, Ry )
% SENSE reconstruction algorithm
% 
% INPUT PARAMETERS:
%  aliasedImages = SIZEXred x SIZEYred x Nc
%                = aliased images at reduced FOV
%  sensProfiles  = SIZEX    x SIZEY    x Nc
%                = maps of the coil sensitivities
%  Rx = acceleration_factor_x
%                = acceleration factor on x
%  Ry = acceleration_factor_y
%                = acceleration factor on y
%  
% OUTPUT PARAMETERS:
%  reconstruction = SIZEX x SIZEY
%                 = SENSE reconstruction


% Get dimensions
[SIZEX,    SIZEY,    Nc] = size(sensProfiles);
[SIZEXred, SIZEYred, Nc] = size(aliasedImages);

% Verify consistency between 
% reduced FOV size and acceleration factor
if (ceil(SIZEX / Rx) ~= SIZEXred) || (ceil(SIZEY / Ry) ~= SIZEYred)
    disp('Aliased Images and Acceleration Factor conflict.');
    return
end

% Reconstruction FOV
reconstruction = zeros(SIZEX, SIZEY);

% Calculate number of replicates
NRX = Rx; %2 * ceil((Rx-1)/2) + 1; % verify why you get 3 
NRY = Ry; %2 * ceil((Ry-1)/2) + 1;

% SENSE reconstruction applied for each voxel in the reduced FOV
for i = 1:SIZEXred
    for j = 1:SIZEYred
   
        % 1. Create aliased signal matrix
        % Every element in the aliased signal matrix
        % will be made up of the signal coming from each
        % coil at the current position
        Saliased = aliasedImages(i, j, :);
        Saliased = Saliased(:);

        % 2. Create sensitivity profile matrix
        % For each replicate (NRX or NRY)
        for nrx = 0 : NRX-1
            for nry = 0 : NRY-1
                
                % Calculate position (y + pL/R) 
                %           y == i, p == nrx
                posx = mod(( (i - 1) + nrx * SIZEXred), SIZEX) + 1;
                posy = mod(( (j - 1) + nry * SIZEYred), SIZEY) + 1;
                
                % Get coils sensitivities for all coils, one position
                % Basically one column vector from the Bcoils matrix
                Bcol = sensProfiles(posx, posy, :);
                Bcol = Bcol(:);
                
                % Every line in the sensitivity profile matrix
                % will be made up of 
                % a coil's sensitivity at each replicate
                if ( (nrx == 0) && (nry == 0))
                    Bcoils = Bcol;
                else
                    Bcoils = [Bcoils Bcol];
                end
                
            end
        end
        
        % 3. Perform reconstruction
        mSignal = pinv(Bcoils) * Saliased;
        
        
        % 4. Populate final reconstruction matrix
        if NRX > 1
            for nrx = 0 : NRX-1
                reconstruction(i+nrx*SIZEXred,j) = mSignal(nrx+1);
            end
        end
        
        if NRY > 1
            for nry = 0 : NRY-1
                reconstruction(i,j+nry*SIZEYred) = mSignal(nrx+nry+2);
            end
        end
        
    end
end



















end

