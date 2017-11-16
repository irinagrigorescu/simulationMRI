% % % % IRINA GRIGORESCU
% % % % This script is to test SENSE 2D

clear all; close all; clc;
addpath ../../helpers

%%
% % An example with a picture
A = imread('~/Pictures/lena.bmp');
[NSIZE, MSIZE, PSIZE] = size(A);
A = A(1:NSIZE,1:NSIZE,1);
A = double(A);

mask = ones(NSIZE,NSIZE);
mask(2:2:end,:) = 0;

% Create Coils
C1 = zeros(NSIZE,NSIZE);
C2 = zeros(NSIZE,NSIZE);
C3 = zeros(NSIZE,NSIZE);
C4 = zeros(NSIZE,NSIZE);
C5 = zeros(NSIZE,NSIZE);
C6 = zeros(NSIZE,NSIZE);
for i = 1:NSIZE
    for j = 1:NSIZE
%         C1(i,j) = exp(i/200)./20;
%         C2(i,j) = exp(-i/200)./20;
%         C3(i,j) = exp(i/100)./20;
%         C4(i,j) = exp(-i/100)./20;
        C1(i,j) = exp(i/350+j/200)./20;
        C2(i,j) = exp(-i/400+j/220)./20;
        C3(i,j) = exp(i/400-j/220)./20;
        C4(i,j) = exp(-i/350-j/200)./20;
    end
end

C1 = abs(C1./max(max(C1)));
C2 = abs(C2./max(max(C2)));
C3 = abs(C3./max(max(C3)));
C4 = abs(C4./max(max(C4)));

sensProfiles = cat(3, C1, C2, C3, C4);
figure, 
subplot(2,2,1), imagesc(C1);
subplot(2,2,2), imagesc(C2);
subplot(2,2,3), imagesc(C3);
subplot(2,2,4), imagesc(C4);
colormap gray

%
% % % Alias images
R = 2;
imAl1 = fft2(C1.*A); imAl1 = ifft2(imAl1(1:R:end,:));
imAl2 = fft2(C2.*A); imAl2 = ifft2(imAl2(1:R:end,:));
imAl3 = fft2(C3.*A); imAl3 = ifft2(imAl3(1:R:end,:));
imAl4 = fft2(C4.*A); imAl4 = ifft2(imAl4(1:R:end,:));
figure, 
subplot(2,2,1), imagesc(imAl1), axis equal;
subplot(2,2,2), imagesc(imAl2), axis equal;
subplot(2,2,3), imagesc(imAl3), axis equal;
subplot(2,2,4), imagesc(imAl4), axis equal;
colormap gray

aliasedIm = cat(3, ...
    imAl1, ...
    imAl2, ...
    imAl3, ...
    imAl4);

% % % Reconstruction
[reconstruction] = sense2D( aliasedIm, sensProfiles, R, 1 );
figure, imagesc(abs(reconstruction)), colormap gray;







