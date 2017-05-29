function Rx = Rotx(theta, flag)
% % % % Rotation about x 
% % % % anti-clockwise fashion

if nargin < 2
    flag = 'aclock';
end

% For clockwise, change theta sign
if strcmp(flag, 'clock') 
    theta = -theta;
end

Rx = [ 1         0            0     ; 
       0   cos(theta)  -sin(theta)  ;
       0   sin(theta)   cos(theta) ];

    
end