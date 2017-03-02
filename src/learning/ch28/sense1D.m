function [ reconstruction ] = sense1D( aliasedImages, sensProfiles, R)
% SENSE reconstruction algorithm
% 
% INPUT PARAMETERS:
%  aliasedImages = SIZEred x Nc
%                = aliased function at reduced FOV
%  sensProfiles  = SIZE    x Nc
%                = maps of the coil sensitivities
%  R             = acceleration factor
%  
% OUTPUT PARAMETERS:
%  reconstruction = SIZE
%                 = SENSE reconstruction


% Get dimensions
[SIZE,    Nc] = size(sensProfiles);
[SIZEred, Nc] = size(aliasedImages);

% Verify consistency between 
% reduced FOV size and acceleration factor
if (ceil(SIZE / R) ~= SIZEred)
    disp('Aliased Images and Acceleration Factor conflict.');
    return
end

% Reconstruction FOV
reconstruction = zeros(SIZE, 1);

% Calculate number of replicates
NR = 2 * ceil((R-1)/2) ;

% SENSE reconstruction applied for each voxel in the reduced FOV
for i = 1:SIZEred
   
        % 1. Create aliased signal matrix
        % Every element in the aliased signal matrix
        % will be made up of the signal coming from each
        % coil at the current position
        Saliased = aliasedImages(i, :);
        Saliased = Saliased(:);

        % 2. Create sensitivity profile matrix
        % For each replicate (NR)
        for nr = 0 : NR-1
                
            % Calculate position (y + pL/R) 
            %           y == i, p == nr
            pos = mod(( (i - 1) + nr * SIZEred), SIZE) + 1;

            % Get coils sensitivities for all coils, one position
            % Basically one column vector from the Bcoils matrix
            Bcol = sensProfiles(pos, :);
            Bcol = Bcol(:);

            % Every line in the sensitivity profile matrix
            % will be made up of 
            % a coil's sensitivity at each replicate
            if (nr == 0)
                Bcoils = Bcol;
            else
                Bcoils = [Bcoils Bcol];
            end
                
        end
        
        % 3. Perform reconstruction
        mSignal = pinv(Bcoils) * Saliased;
        
        % 4. Populate final reconstruction matrix
        reconstruction(i) = mSignal(1,1);
        reconstruction(i+SIZEred) = mSignal(2,1);
        
    
end



















end

