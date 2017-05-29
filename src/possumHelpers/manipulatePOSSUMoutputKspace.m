% % % % IRINA GRIGORESCU
% % % % This script will manipulate the output of a possum simulation


% Prerequisites
addpath ../helpers/
addpath ~/Tools/MatlabStuff/matlabnifti/

%% Load files
filename = '/Users/irina/simdir/outKspace/kspace_';

TRs = 25;

% Allocating memory
kspaceData = struct;
kspaceData.real = zeros(64,64,TRs);
kspaceData.imag = zeros(64,64,TRs);
kspaceData.sumreal = zeros(1,TRs);
kspaceData.sumimag = zeros(1,TRs);

for tr = 1:TRs
    % Load images
    kspaceDataReal = load_nii(...
            [filename, num2str(tr), '_real.nii.gz'], ...
              [], [], [], [], [], 0.5);
	
	kspaceDataImag = load_nii(...
            [filename, num2str(tr), '_imag.nii.gz'], ...
              [], [], [], [], [], 0.5);
          
	% Save data in variables
    kspaceData.real(:,:,tr) = kspaceDataReal.img;
    kspaceData.imag(:,:,tr) = kspaceDataImag.img;
    kspaceData.sumreal(1,tr) = sum(sum(kspaceDataReal.img));
    kspaceData.sumimag(1,tr) = sum(sum(kspaceDataImag.img));
end

%% Plot their absolute value
% just middle kspace
% posx = 33; posy = 33;



signal = abs (    kspaceData.sumreal(1,:).^2 + ...
              1i.*kspaceData.sumreal(1,:).^2);

figure
plot(1:TRs, signal)
xlabel('TR index')
ylabel('Signal')

