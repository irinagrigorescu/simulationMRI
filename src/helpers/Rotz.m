function Rz = Rotz(theta, flag)
% % % % IRINA GRIGORESCU
% % % % Rotation about z
% % % % anti-clockwise fashion

if nargin < 2
    flag = 'aclock';
end

% For clockwise, change theta sign
if strcmp(flag, 'clock') 
    theta = -theta;
end

Rz = [ cos(theta) -sin(theta) 0 ; 
       sin(theta)  cos(theta) 0 ;
          0            0      1 ];  
    
end