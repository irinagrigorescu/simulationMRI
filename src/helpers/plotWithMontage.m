function plotWithMontage(tmp,m,n,p,size)
% INPUT:
%   tmp     = 3D volume
%   [m,n,p] = size(tmp)
%   size    = number of columns for window
% OUTPUT:
%   figure

    figure
    for i = 1:p
        tmp2(1:m, 1:n, 1, i) = tmp(1:m,1:n,i)./max(tmp(:));
    end
    montage(tmp2, 'Size', [NaN size]), colormap jet
    clear tmp
    clear tmp2
end