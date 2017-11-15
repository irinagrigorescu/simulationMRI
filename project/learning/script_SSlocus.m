% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % Irina Grigorescu
% % % % Date created: 10-08-2017
% % % % Date updated: 10-08-2017
% % % % 
% % % % Based on several publications such as:
% % % % Scheffler, On the transient phase of balanced SSFP sequences
% % % % Freeman, Phase and intensity anomalies in fourier transform NMR
% % % % Zur, An analysis of fast imaging sequences with SS transverse magnetization refocusing
% % % % Dharmakumar, Understanding steady-state free precession: A geometric perspective
% % % % 
% % % % I have come to the conclusion that the SS locus of the
% % % % magnetization is interesting to investigate in terms of it's 
% % % % geometric representation.
% % % % This script is my attempt at plotting and understanding its shape.

addpath ../helpers

% % % % Prerequisites for simulation:
T1    = 950;
T2    = T1;%100;
TR    =  20;
theta = deg2rad( 45);
beta  = deg2rad(-180:180);
E1    = exp(-TR/T1);
E2    = exp(-TR/T2);
M0    =   1;

% % % % The equilibrium solution for an SSFP signal (meaning that the
% resonance offset angle is taken into account) comes from:
% 
% 1. Writing up the component-wise magnetization vectors during a TR block
%    as equations depending on the same components immediately after the RF
%    pulse: Mx(n,t') = [ Mx+(n) cos b(t')  + My+(n) sin b(t')] e^(-t'/T2)
%           My(n,t') = [-Mx+(n) sin b(t')  + My+(n) cos b(t')] e^(-t'/T2)
%           Mz(n,t') =   Mz+(n) e^(-t'/T1) + M0 (1 - e^(-t'/T1))
%    where t' = t - nTR (relative time during block)
% 2. Write M+ as M+(n  ) = Rx(theta) M-(n)
% 3. And   M- as M-(n+1) = D(TR) M+(n)   + M0(1-E1)z
% 4. And         M-(n  ) = D(TR) M+(n-1) + M0(1-E1)z
% 5. And set     M-(n  ) = M-(n+1)
% You get the solutions for the M+(inf) and M-(inf)
% where d = 
d = (1 - E1 * cos(theta)) .* (1 - E2 .* cos(beta)) - ...
   E2 * (E1 - cos(theta)) .* (E2 - cos(beta));

% Component-wise, these solutions are:
MpxSS = M0*(1-E1) .* (E2*sin(theta).*sin(beta)) ./ d;
MpySS = M0*(1-E1) .* (sin(theta)*(1-E2.*cos(beta))) ./ d;
MpzSS = M0*(1-E1) .* (E2*(E2 - cos(beta)) + (1 - E2.*cos(beta))*cos(theta)) ./ d;

% Plot this structure
figure
createAxis(1);
view(127,17);
scatter3(MpxSS, MpySS, MpzSS, '.');
axis square

% To describe this structure in terms of an ellipse we can write the
% ellipse equation: (x-x0)^2/a^2 + (y-y0)^2/b^2 = 1 where:
%                   (x0, y0) are the ellipse's center and
%                   ( a,  b) are the ellipse's 'radii'
% As in our example, the RF pulse tips the magnetization about the x-axis,
% it follows that x0 = 0 for the coordinates of the MpxSS will be centered
% around 0.
% If we describe the ellipse as having La and Lb as the major and minor
% diameters of the ellipse we can write:
% x^2/(La/2)^2 + (y-(Lb/2))^2/(Lb/2)^2 = 1

% Further analysis, such as setting d Mx / d beta = 0 and d My / d beta = 0
% will yield the extreme values for the locus in both x and y directions.
% which will give La = Lb for T1 = T2










