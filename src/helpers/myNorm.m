function [ M_norm ] = myNorm( M )
% % % % IRINA GRIGORESCU 01/11/2016
% % % % Normalizes an array

    max_val = max(abs(M(:)));
    
    if max_val ~= 0
        M_norm = M./max_val;
    else
        M_norm = M;
    end
    
end

