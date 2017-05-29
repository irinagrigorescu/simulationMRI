function Ry = Roty(theta, flag)
% % % % Rotation about y
% % % % anti-clockwise fashion

if nargin < 2
    flag = 'aclock';
end

% For clockwise, change theta sign
if strcmp(flag, 'clock') 
    theta = -theta;
end

Ry = [ cos(theta)  0   sin(theta) ; 
          0        1         0    ;
      -sin(theta)  0   cos(theta) ];
    
end