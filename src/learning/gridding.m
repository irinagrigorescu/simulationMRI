
% Interpolate a 2D scattered data set over a uniform grid
xy = -2.5 + 5*gallery('uniformdata',[200 2],0);
x = xy(:,1); y = xy(:,2);
% x = linspace(-2,2,200)'; 
% y = linspace(-2,2,200)';
v = x.^2 + y.^2 + 25;
[xq,yq] = meshgrid(-2:.2:2, -2:.2:2);
vq = griddata(x,y,v,xq,yq);
mesh(xq,yq,vq), hold on, plot3(x,y,v,'o'), hold off


