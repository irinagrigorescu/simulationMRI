function deletePlots(plotToDelete, varargin)
% % This function deletes all plots sent as params
% % IRINA GRIGORESCU

    % at least one parameter is needed and therefore deleted here
    delete(plotToDelete);

    % the rest are deleted here
    nVarargs = length(varargin);
    for k = 1:nVarargs
        delete(varargin{k});
    end
end