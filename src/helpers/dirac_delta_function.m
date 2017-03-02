function [ax] = dirac_delta_function(x, ths)
% % % % IRINA GRIGORESCU 
% % % % 25/10/2016
% % % % The Dirac delta function is +inf when x ~= 0
% % % % 
% % % % In our discrete case we return 1 when the argument is smaller than
% a threshold as a common use for Dirac Delta is d(x-a), where it would be
% +infinity when x-a == 0, but here we do x-a < epsilon.

    ax = abs(x) < ths;    
    
end
