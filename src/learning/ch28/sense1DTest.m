% % % % IRINA GRIGORESCU
% % % % This script is to test SENSE 1D

clear all; close all; clc;
addpath ../../helpers

%% 1. Initial 1D signal
% Set some parameters (hard coded for this example)
Ny = 64; % Number of points in the signal
L  = 12;  % field-of-view
Lmin = -L/2; % negative side of the plot 
Lmax =  L/2; % positive side of the plot

% Create the initial signal
y  = linspace(-L/2, L/2, Ny); % y spans from -L/2 to L/2
My = zeros(1, Ny);      % My initialised as 0 for each of the Ny points

% M(y) will take these values 
for i = 1:Ny
    
    % 1. 0.5 for -5 < y < 0
    if y(i) > Lmin+1 && y(i) < 0
        My(i) = 0.5;
    end
    
    % 2. 1.0 for  0 < y < 5
    if y(i) > 0 && y(i) < Lmax-1
        My(i) = 1.0;
    end
    
end

% Plotting the signal
figure, 

% 1. Plot the initial My signal
subplot(2,1,1)
colorsPlot = {'--d','-o','--*'};
% % % Repeat the signal to simulate behaviour of Fourier periodicity
for p = -1:1
    % Periodicity calculated as y+p*L,
    % where p is here hardcoded but can be calculated as:
    %
    % p = -N(A/L) to N(A/L) (see Ch28, eq 28.8), where
    %
    % N(A/L) = the least non-negative integer larger than or equal to
    % (A/L-1)/2
    %
    % Obs: Here, we make the assumption that A = L
    plot(y+p*L, My, [colorsPlot{p+2}, 'b']); hold on
end
title(['Initial 1D signal repeated every L = ', num2str(L), ' steps']);
% Plotting labels for the x axis
ax = gca;
axis([-L-L/2-1, L+L/2+1, -.5, 1.2]), grid on
ax.XTick = [-3*L/2, -L, -L/2, 0, L/2, L, 3*L/2];
ax.XTickLabel = {'-3L/2', '-L', '-L/2', '0', 'L/2', 'L', '3L/2'};
xlabel('y (s)'); ylabel('M(y)');

% 2. Fourier transform the signal and back
%    and plot it
MyTrans = ifft(fft(My));
subplot(2,1,2)
colorsPlot = {'--d','-o','--*'};
% % % Repeat the signal to simulate behaviour of Fourier periodicity
for p = -1:1
    % Periodicity calculated as y+p*L,
    % where p is here hardcoded but can be calculated as:
    %
    % p = -N(A/L) to N(A/L) (see Ch28, eq 28.8), where
    %
    % N(A/L) = the least non-negative integer larger than or equal to
    % (A/L-1)/2
    %
    % Obs: Here, we make the assumption that A = L
    plot(y+p*L, real(MyTrans), [colorsPlot{p+2}, 'b']); hold on
end
title(['Initial (Fourier transformed) 1D signal repeated every L = ', num2str(L), ' steps']);
% Plotting labels for the x axis
ax = gca;
axis([-L-L/2-1, L+L/2+1, -.5, 1.2]), grid on
ax.XTick = [-3*L/2, -L, -L/2, 0, L/2, L, 3*L/2];
ax.XTickLabel = {'-3L/2', '-L', '-L/2', '0', 'L/2', 'L', '3L/2'};
xlabel('y (s)'); ylabel('iFFT(FFT(M(y)))');


%% 2. Reduce FOV from L to L/R
R = 2; % acceleration factor

% 1. Plot the initial My signal
figure
subplot(2,1,1)
colorsPlot = {'--d','-o','--*'};
% % % Repeat the signal to simulate behaviour of Fourier periodicity
for p = -1:1
    % Periodicity calculated as y+p*L,
    % where p is here hardcoded but can be calculated as:
    %
    % p = -N(A/L) to N(A/L) (see Ch28, eq 28.8), where
    %
    % N(A/L) = the least non-negative integer larger than or equal to
    % (A/L-1)/2
    %
    % Obs: Here, we make the assumption that A = L
    plot(y+p*L, My, [colorsPlot{p+2}, 'b']); hold on
end
title(['Initial 1D signal repeated every L = ', num2str(L), ' steps']);
% Plotting labels for the x axis
ax = gca;
axis([-L-L/2-1, L+L/2+1, -.5, 1.2]), grid on
ax.XTick = [-3*L/2, -L, -L/2, 0, L/2, L, 3*L/2];
ax.XTickLabel = {'-3L/2', '-L', '-L/2', '0', 'L/2', 'L', '3L/2'};
xlabel('y (s)'); ylabel('M(y)');

% 2. Plot the reduced FOV signal
% The reduced FOV signal will show aliasing which can be viewed as a 
% juxtaposition of signals as the reduced FOV brings closer together the 
% Fourier replicates
subplot(2,1,2)
colorsPlot = {'--dr','-ob','--*g'};
% % % Repeat the signal to simulate behaviour of Fourier periodicity
for p = -1:1
    % Periodicity calculated as y+p*L,
    % where p is here hardcoded but can be calculated as:
    %
    % p = -N(AR/L) to N(AR/L) (see Ch28, eq 28.8), where
    %
    % N(AR/L) = the least non-negative integer larger than or equal to
    % (AR/L-1)/2
    %
    % Obs: Here, we make the assumption that A = L
    plot(y+p*L/R, My, [colorsPlot{p+2}, '']); hold on
end
title(['Initial 1D signal repeated every L/R = ', ...
    num2str(L), ' / ', num2str(R), ' = ', num2str(L/R), ' steps']);
% Plotting labels for the x axis
ax = gca;
axis([-L-L/2-1, L+L/2+1, -.5, 1.2]), grid on
ax.XTick = [-3*L/2, -L, -L/2, 0, L/2, L, 3*L/2];
ax.XTickLabel = {'-3L/2', '-L', '-L/2', '0', 'L/2', 'L', '3L/2'};
xlabel('y (s)'); ylabel('M(y)');

%% 3. Simulate aliasing with k-space undersampling
%    Aliasing occurs when the reduced FOV is smaller than the object's size
%    This can be simulated as such:
%    3.1 The k-space step is increased with the R factor making it 
%        dKnew = R * dK
%        Because dK = 1/FOV => dKnew = 1 / (FOV/R) => FOVnew = FOV/R
%    3.2 The k-space range is kept the same, meaning that:
%        dx and dy will be kept the same size
%        dx = 1/(2 kxmax) and dy = 1/(2 kymax)
%        Obs: This is a bit evident as dx = 1 / (dkx Nx) and
%             we increase dkx but decrease Nx
MyTrans = fft(My);
yred = linspace(-L/(R*2), L/(R*2), Ny/R); % y spans from -L/(R*2) to L/(R*2) (image space)


figure,
% Plot the initial M(y) signal
subplot(3,1,1)
colorsPlot = {'--d','-o','--*'};
% % % Repeat the signal to simulate behaviour of Fourier periodicity
for p = 0:0
    % See above for explanation
    %
    % Obs: Here, we make the assumption that A = L and we chose to plot
    % only one repetition of the signal
    plot(y+p*L, My, [colorsPlot{p+2}, 'b']); hold on
end
title(['Initial 1D signal repeated every L = ', num2str(L), ' steps']);
% Plotting labels for the x axis
ax = gca;
axis([-L-L/2-1, L+L/2+1, -.5, 1.2]), grid on
ax.XTick = [-3*L/2, -L, -L/2, 0, L/2, L, 3*L/2];
ax.XTickLabel = {'-3L/2', '-L', '-L/2', '0', 'L/2', 'L', '3L/2'};
xlabel('y (s)'); ylabel('M(y)');

% Plot k-space and undersampled k-space
subplot(3,1,2)
colorsPlot = {'--d','-o','--*'};
% % % Repeat the signal to simulate behaviour of Fourier periodicity
for p = 0:0
    % See above for explanation
    %
    % Obs: Here, we make the assumption that A = L and we chose to plot
    % only one repetition of the signal
    plot(y+p*L, abs(fftshift(MyTrans)), [colorsPlot{p+2}, 'b']); hold on
    plot(y(1:R:end)+p*L, abs(fftshift(MyTrans(1:R:end))), [colorsPlot{3}, 'r']); hold on
end
title('Signal in K-space');
% Plotting labels for the x axis
ax = gca;
axis([-L-L/2-1, L+L/2+1, min(MyTrans)-1, max(MyTrans)+1]), grid on
ax.XTick = [-L/2, 0, L/2];
ax.XTickLabel = {'-Kmax', '0', '+Kmax'};
xlabel('k (Hz)'); ylabel('abs (FFT( M(y) ))');
legend('R = 1', 'R = 2');


% Plot the inverse fourier transform
subplot(3,1,3)
colorsPlot = {'--d','-o','--*'};
% % % Repeat the signal to simulate behaviour of Fourier periodicity
for p = -1:1
    % See above for explanation
    %
    % Obs: Here, we make the assumption that A = L and we chose to plot
    % only one repetition of the signal
    %plot(y+p*L, abs(ifft(MyTrans)), [colorsPlot{p+2}, 'b']); hold on
    plot(y+p*L/R, abs(ifft(MyTrans)), [colorsPlot{p+2}, 'b']); hold on
end
for p = 0:1
    plot(yred-L/(2*R)+p*L/R, abs(ifft(MyTrans(1:R:end))), [colorsPlot{3}, 'r']); hold on
end
title('Fourier transform of k-space');
% Plotting labels for the x axis
ax = gca;
axis([-L-L/2-1, L+L/2+1, -.5, 2]), grid on
ax.XTick = [-3*L/2, -L, -L/2, -L/(2*R), 0, L/(2*R), L/2, L, 3*L/2];
ax.XTickLabel = {'-3L/2', '-L', '-L/2', '-L/(2R)', '0', 'L/(2R)', 'L/2', 'L', '3L/2'};
xlabel('y (s)'); ylabel('abs(iFFT(FFT(M(y))))');
legend('R = 1', 'R = 1', 'R = 1', 'R = 2');


%% This shows the non-symmetry of k-space 
%  when subsampling is not done with acceleration factors powers of 2.
Rmax = 8;
figure
for i = 1:Rmax
    ll = abs(fftshift(MyTrans(1:i:end)));
    
    subplot(Rmax,2,(2*i)-1); 
    imagesc(ll); title(['R = ', num2str(i)]);
    ax = gca;
    ax.XTickLabel = {};
    ax.YTickLabel = {};
    
    subplot(Rmax,2,2*i); 
    plot(ll,'*--'); 
    ax = gca;
    ax.XTickLabel = {};
    ax.YTickLabel = {};
end


return
%% 2. Coils profiles
c1 = zeros(1, N);
c2 = zeros(1, N);
for i = 1:N
    c2(i) = c2(i) + i^(.8);
    c1(N-i+1) = c2(i)+0.02;
end
c = [c1 ; c2]';
c = c./max(max(c));

% Plotting the coils's profiles
for p = -1:1
    plot(Lmin+p*L : Lmax+p*L, c(:,1), 'r'); hold on
    plot(Lmin+p*L : Lmax+p*L, c(:,2), 'g'); hold on
end
title('1D signal and coils');
legend('si','si','si','c1','c2');

axis([Lmin-L-1, Lmax+L+1, -.5, 1.2]), grid on


%% 3. Signal coming from each coil
sc1 = si.*c(:,1);
sc2 = si.*c(:,2);

subplot(2,1,2)
for p = -1:1
    plot(Lmin+p*L : Lmax+p*L, sc1, [colorsPlot{p+2},'r']); hold on
    plot(Lmin+p*L : Lmax+p*L, sc2, [colorsPlot{p+2},'g']); hold on
end
title('Unaliased Signal coming from each coil');
legend('Sc1','Sc2');

ax = gca;
axis([Lmin-L-1, Lmax+L+1, -.5, 1.2]), grid on
ax.XTick = [1-L, 1, Lmax, Lmax+L];
ax.XTickLabel = {'1-L', '1', 'L+1', '(L+1)+L'};

%% 4. Alias the signal
% Let's alias the signal
sac1 = fft(sc1); sac1red = ifft(sac1(1:R:end));
sac2 = fft(sc2); sac2red = ifft(sac2(1:R:end));
sacred = [sac1red sac2red];

% Plot unaliased and aliased signals
figure,
subplot(3,1,1)
for p = -1:1
    plot(Lmin+p*L : Lmax+p*L, sc1, [colorsPlot{p+2},'r']); hold on
    plot(Lmin+p*L : Lmax+p*L, sc2, [colorsPlot{p+2},'g']); hold on
end
title('Unaliased Signal coming from each coil');
legend('Sc1','Sc2');

ax = gca;
axis([Lmin-L-1, Lmax+L+1, -.5, 1.2]), grid on
ax.XTick = [1-L, 1, Lmax, Lmax+L];
ax.XTickLabel = {'1-L', '1', 'L+1', '(L+1)+L'};

% Plot the juxstaposed signals
subplot(3,1,2)
for p = -1:1
    plot(Lmin+p*L/R : Lmax+p*L/R, sc1, 'r'); hold on
    plot(Lmin+p*L/R : Lmax+p*L/R, sc2, 'g'); hold on
end
title('Reduced FOV causes aliasing');
legend('Sc1','Sc2');

ax = gca;
axis([Lmin-L-1, Lmax+L+1, -.5, 1.2]), grid on
ax.XTick = [1-L/1, 1-L/R, 1, Lmin+L/R, Lmin+L, Lmin+L+L/R, Lmin+L+L];
ax.XTickLabel = {'1-L','1-L/R', '1', '1+L/R', '1+L', '(L+1)+L/R', '(L+1)+L'};

% Plot the aliased signals coming from each coil
subplot(3,1,3)
Lminn = Lmin;%(Lmin+L/R)-(L/R)/2;
Lmaxn = Lmin+L/R; %(Lmax-L/R)+(L/R)/2;
for p = 0:0
   plot(Lminn+p*L/R : Lmaxn+p*L/R , sac1red, 'r'); hold on
   plot(Lminn+p*L/R : Lmaxn+p*L/R , sac2red, 'g'); hold on
end
title('Aliased Signal coming from each coil');
legend('Sac1','Sac2');

ax = gca;
axis([Lmin-L-1, Lmax+L+1, -.5, 2.2]), grid on
ax.XTick = [1-L/1, 1-L/R, 1, Lmin+L/R, Lmin+L, Lmin+L+L/R, Lmin+L+L];
ax.XTickLabel = {'1-L','1-L/R', '1', '1+L/R', '1+L', '(L+1)+L/R', '(L+1)+L'};


return
%%
% Reconstruction
reconstructed1D = sense1D( sacred, c, 2);
subplot(2,2,4)
plot(lmin:lmax, abs(reconstructed1D)); 
title('Reconstruction');
legend('Srec');
axis([lmin-1,lmax+1,-.5,1.2]), grid on



%% Some other visualisations

% How aliasing occurs
figure
subplot(2,1,1)
plot(lmin         : lmax           , sc1 ,  '-b'); hold on
plot(lmin - L/2   : 0              , sc1 ,  '-r'); hold on
plot(  0          : lmax + L/2 - 1 , sc1 , '--g'); hold on
%plot(-5:1, real(sacred(:,1))); hold on
%plot( 1:7, real(sacred(:,1))); hold on
title('Aliasing for the signal coming from coil 1');
legend('Sc1');
axis([-N-1, N+1,-.5,1.2]), grid on
ax = gca;
ax.XTick = [-N+1, -N/2, 0, N/2, N-1];
ax.XTickLabel = {'-L', '-L/2', '0', 'L/2', 'L'};

subplot(2,1,2)
plot(lmin:lmax, [ real(sacred(4:end,1)) ; real(sacred(:,1)) ; real(sacred(1:3,1)) ] ); hold on
title('Aliased Signal coming from coil 1');
axis([-N-1, N+1,-.5,1.2]), grid on
ax = gca;
ax.XTick = [-N+1, -N/2, 0, N/2, N-1];
ax.XTickLabel = {'-L', '-L/2', '0', 'L/2', 'L'};














