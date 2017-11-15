function func_saveToVideo(videoClass, figHandle, videoFlag, toClose)
% 
% This function can be used to save your animation into a video file
% 
% Input:
%   videoClass = video class object
%   figHandle = the handle to the figure for saving
%   videoFlag = struct with 2 entries
%               .flag = 0 (don't save) | 1 (save)
%               .nameVideo = 'name_of_video_file' (e.g. 'spinPrecession')
%   toClose   = 1 (will close video handle) 
%               0 (will not close video handle)
% 
% Author: Irina Grigorescu, irina.grigorescu.15@ucl.ac.uk
%                           irinagry@gmail.com
% 

% To turn on/off your save to video file functionality
% I have chosen to add this flag which can be easily set to 0
% from outisde of the function
% Bare in mind that saving to a video file will make your simulation run
% slower

if videoFlag.flag == 1
    
    % Get current frame from your figure handle
    frame = getframe(figHandle);               % #need this for video
    % Write the frame to the video object
    writeVideo(videoClass, frame)   
    
end

% Finally, close the video handle
if toClose == 1
    close(videoClass)                      % #need this for video
end

end