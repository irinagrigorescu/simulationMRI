function plotWithMontage(tmp,m,n,p,size)
    figure
    for i = 1:p
        tmp2(1:m, 1:n, 1, i) = tmp(1:m,1:n,i)./max(tmp(:));
    end
    montage(tmp2, 'Size', [NaN size]), colormap jet
    clear tmp
    clear tmp2
end