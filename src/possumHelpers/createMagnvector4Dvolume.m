% % % % IRINA GRIGORESCU
% % % % Creating a magnetisation vector xyz components 4D volume


% Prerequisites
addpath ../helpers/
addpath ~/Tools/MatlabStuff/matlabnifti/

%% Load file
filename = '~/simdir/brain.nii.gz';
iniData = load_nii(filename, ...
              [], [], [], [], [], 0.5);

sizeOfData = iniData.original.hdr.dime.dim(2:5);

%% See brain
% brainTissues = {'GM', 'WM', 'CSF'};
% for i = 1:sizeOfData(4)
%     plotWithMontage(squeeze(outputData.img(:,:,:, i)), ...
%                     sizeOfData(1), ...
%                     sizeOfData(2), ...
%                     sizeOfData(3), ...
%                     15);
% 	title(brainTissues{i});
%     pause(5)
% end

%% Find a place in the brain where it's only WM
% for i = 1:sizeOfData(1)
%     for j = 1:sizeOfData(2)
%         for k = 1:sizeOfData(3)
%             if iniData.img(i,j,k, 2) == 1 && ...
%                 (iniData.img(i,j,k, 1) == 0 && iniData.img(i,j,k, 3) == 0)
%                 i 
%                 j
%                 k
%                 iniData.img(i,j,k,:)
%                 return
%             end
%         end
%     end
% end

%% Create empty 4D volume of the same size and properties as the 
%  loaded file
outputData = iniData;

% Increase 4th component to a value
fourthDim = 3; % the x/y/z components
outputData.original.hdr.dime.dim(5) = fourthDim;
outputData.hdr.dime.dim(5) = fourthDim;
newSizeOfData = sizeOfData;
newSizeOfData(4) = fourthDim;
outputData.img = zeros(newSizeOfData);
% outputData.img( 91, 109, 49, 3) = 1;
% outputData.img( 73,  31, 47, 3) = 1;
outputData.img( 23, 95, 75, 3) = -1; % here it's only WM
squeeze(iniData.img( 23, 95, 75, :))
% outputData.img(100:110, 30:40, 47, 3) = 1;
outputData.fileprefix = '~/simdir/brainmagn_1.nii.gz';
save_nii(outputData, '~/simdir/brainmagn_1.nii.gz');











