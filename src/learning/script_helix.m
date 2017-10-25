% Irina Grigorescu
% This script looks at how a helix can be represented in terms of
% Euler formula
% 



k = 100;

figure
for r = 0:0.001:1
    scatter3(real(exp(-1i.*k.*r)), imag(exp(-1i.*k.*r)), r, 'b.'); 
    hold on
    xlim([-1 1]); ylim([-1 1]); zlim([0 1]);
    %pause(0.01)
    drawnow
end






