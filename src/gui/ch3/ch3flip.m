function varargout = ch3flip(varargin)
% CH3FLIP MATLAB code for ch3flip.fig
%      CH3FLIP, by itself, creates a new CH3FLIP or raises the existing
%      singleton*.
%
%      H = CH3FLIP returns the handle to a new CH3FLIP or the handle to
%      the existing singleton*.
%
%      CH3FLIP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CH3FLIP.M with the given input arguments.
%
%      CH3FLIP('Property','Value',...) creates a new CH3FLIP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ch3flip_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ch3flip_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ch3flip

% Last Modified by GUIDE v2.5 28-Nov-2016 20:42:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ch3flip_OpeningFcn, ...
                   'gui_OutputFcn',  @ch3flip_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% % % % IRINA GRIGORESCU
% % % % DATE: 28-Nov-2016
% % % % GUI for Chapter 3 stuff

% --- Executes just before ch3flip is made visible.
function ch3flip_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ch3flip (see VARARGIN)

% Choose default command line output for ch3flip
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ch3flip wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% % % % Add path to algorithms
addpath ../../algs/ch3/
addpath ../../helpers/


% --- Outputs from this function are returned to the command line.
function varargout = ch3flip_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% % % % % % HANDLES
% handles = 
% 
%              figure1: [1x1 Figure]
%             uipanel6: [1x1 Panel]
%       uibuttongroup2: [1x1 ButtonGroup]
%       uibuttongroup1: [1x1 ButtonGroup]
%               runSim: [1x1 UIControl]
%             uipanel3: [1x1 Panel]
%             uipanel4: [1x1 Panel]
%             uipanel1: [1x1 Panel]
%                axes1: [1x1 Axes]
%           speedOfSim: [1x1 UIControl]
%           rotateB1mx: [1x1 UIControl]
%           rotateB1my: [1x1 UIControl]
%            rotateB1x: [1x1 UIControl]
%            rotateB1y: [1x1 UIControl]
%             rotframe: [1x1 UIControl]
%             labframe: [1x1 UIControl]
%              B0label: [1x1 UIControl]
%                B0val: [1x1 UIControl]
%           updateGyro: [1x1 UIControl]
%           gammalabel: [1x1 UIControl]
%             gammaval: [1x1 UIControl]
%        updateRFangle: [1x1 UIControl]
%         updateRFtime: [1x1 UIControl]
%             updateB1: [1x1 UIControl]
%     RFflipanglelabel: [1x1 UIControl]
%       RFflipangleval: [1x1 UIControl]
%          RFtimelabel: [1x1 UIControl]
%            RFtimeval: [1x1 UIControl]
%              B1label: [1x1 UIControl]
%                B1val: [1x1 UIControl]
%               output: [1x1 Figure]
% % % % % % END HANDLES

% --- Executes on button press in runSim.
function runSim_Callback(hObject, eventdata, handles)
% hObject    handle to runSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles
% % % % Create figure nicely
handles.figure1;
cla reset;
% Create static coordinate system
[xax, yax, zax] = createAxis();
% Viewer position
view([45 45 10]);

% % % % % Get values given by user
% make it back to real value
gammaValue       = str2double(get(handles.gammaval, 'String')) * (10^8);
% in tesla
b0Value          = str2double(get(handles.B0val, 'String'));
% make it from micro tesla to tesla
b1Value          = str2double(get(handles.B1val, 'String')) * (10^-6);
% make it from ms to sec
rfpulsetimeValue = str2double(get(handles.RFtimeval, 'String')) * (10^-3);
% make it from degrees to radians
%rfflipangleValue = str2double(get(handles.RFflipangleval, 'String')) * pi / 180;
% get speed of simulation
speedUpSim = str2double(get(handles.speedOfSim, 'String'));
if speedUpSim < 0
    speedUpSim = 1;
else
    speedUpSim = floor(speedUpSim);
end
set(handles.speedOfSim, 'String', num2str(speedUpSim)); 

% % % % See if the user wants to run the simulation with rotation or
% without
labFrameValue = get(handles.labframe, 'Value');
% % % % See if the user wants B1 to be static along x-axis or y-axis
xplus        = get(handles.rotateB1x,  'Value');
xminus       = get(handles.rotateB1mx, 'Value');
yplus        = get(handles.rotateB1y,  'Value');

rotateB1Axis = (xplus == 1) || (xminus == 1);

switch rotateB1Axis
    % If it is static along x-axis
    case 1
        if xplus == 1
            b1axis = '+x'; 
        else
            b1axis = '-x';
        end
    % If it is static along y-axis
    case 0
        if yplus == 1
            b1axis = '+y';
        else
            b1axis = '-y';
        end
end
            

% % % % Call appropriate function
switch labFrameValue
    % % % The case where you see all rotations
    case 1
        % Calculate spin movement history
        [ mu, muz, muxy, b1field ] = ...
            flipRotate( gammaValue, b0Value, b1Value, ...
            rfpulsetimeValue, 'no', b1axis );
        
    % % % The case where you are in the rotating frame
    case 0
        % Calculate spin movement history
        [ mu, muz, muxy, b1field ] = ...
            flipRotate( gammaValue, b0Value, b1Value, ...
            rfpulsetimeValue, 'yes', b1axis );
end

% % % % DRAW
N = length(mu);
        
% Plot animation
for i = 2:speedUpSim:N
    % Vector from (0,0,0) to (x,y,z) of updated position
    muvec = [[0 0 0]'  mu(:,i)]; 
    muxyvec = [[0 0 0]'  muxy(:,i)]; 
    muzvec = [[0 0 0]'  muz(:,i)]; 
    b1vec = [[0 0 0]'  b1field(:,i)]; 

    % Plot vector and arrow at the end
    % for mu vector
    muvecplot = plot3(muvec(1,:), muvec(2,:), muvec(3,:), 'r'); hold on
    mutipplot = scatter3(mu(1,i), mu(2,i), mu(3,i), 'r^'); hold on
    % Draw history of tip
    mutiphistplot = scatter3(mu(1,i), mu(2,i), mu(3,i), 'k.'); hold on
    % for xy projection of mu vector
    muxyplot = plot3(muxyvec(1,:), muxyvec(2,:), muxyvec(3,:), '--r'); hold on
    muxytipplot = scatter3(muxy(1,i), muxy(2,i), muxy(3,i), 'r*'); hold on
    % for z projection of mu vector
    muzplot = plot3(muzvec(1,:), muzvec(2,:), muzvec(3,:), '--r'); hold on
    muztipplot = scatter3(muz(1,i), muz(2,i), muz(3,i), 'r*'); hold on
    % for B1 field
    b1vecplot = plot3(b1vec(1,:), b1vec(2,:), b1vec(3,:), '--g'); hold on
    b1tipplot = scatter3(b1field(1,i), b1field(2,i), b1field(3,i), 'ko'); hold on

    % Draw and pause
    drawnow
    pause(0.00001)

    % Delete vector and arrow so that each update looks new
    % Don't delete last one
    if abs(i - N) >= speedUpSim
        deletePlots(muvecplot, mutipplot, ...
                    muxyplot, muxytipplot, ...
                    muzplot, muztipplot, ...
                    b1vecplot, b1tipplot); 
        if labFrameValue == 1
            deletePlots(mutiphistplot);
        end
    end

    % Plot history
    scatter3(b1field(1,i), b1field(2,i), b1field(3,i), 'k.'); hold on
end

% --- END FUNCTION FOR EXECUTION OF RUN

function B0val_Callback(hObject, eventdata, handles)
% hObject    handle to B0val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B0val as text
%        str2double(get(hObject,'String')) returns contents of B0val as a double


% --- Executes during object creation, after setting all properties.
function B0val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B0val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gammaval_Callback(hObject, eventdata, handles)
% hObject    handle to gammaval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gammaval as text
%        str2double(get(hObject,'String')) returns contents of gammaval as a double


% --- Executes during object creation, after setting all properties.
function gammaval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gammaval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B1val_Callback(hObject, eventdata, handles)
% hObject    handle to B1val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B1val as text
%        str2double(get(hObject,'String')) returns contents of B1val as a double


% --- Executes during object creation, after setting all properties.
function B1val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B1val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function RFtimeval_Callback(hObject, eventdata, handles)
% hObject    handle to RFtimeval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RFtimeval as text
%        str2double(get(hObject,'String')) returns contents of RFtimeval as 
%        a double


% --- Executes during object creation, after setting all properties.
function RFtimeval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RFtimeval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RFflipangleval_Callback(hObject, eventdata, handles)
% hObject    handle to RFflipangleval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RFflipangleval as text
%        str2double(get(hObject,'String')) returns contents of RFflipangleval as a double


% --- Executes during object creation, after setting all properties.
function RFflipangleval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RFflipangleval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in updateGyro.
function updateGyro_Callback(hObject, eventdata, handles)
% hObject    handle to updateGyro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % % % When gyro is updated, 
% % % %    keep constant: flip angle, rf pulse time
% % % %    change: B1
% gamma value make it back to 10^6
gammaValue = str2double(get(handles.gammaval, 'String')) * (10^8);
% make it from ms to sec
tauValue   = str2double(get(handles.RFtimeval, 'String')) * (10^-3);
% make it from degrees to radians
thetaValue = str2double(get(handles.RFflipangleval, 'String')) * pi / 180;

% % % % New calculated value
newB1Value = (thetaValue / (tauValue * gammaValue) );
% % % % Change old B1 value with new B1 value
% % % %     multiplying with 10^6 to make it micro
set(handles.B1val,'String',num2str(newB1Value * (10^6)));




% --- Executes on button press in updateB1.
function updateB1_Callback(hObject, eventdata, handles)
% hObject    handle to updateB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % % % When B1 is changed
% % % %     keep constant: RF pulse time
% % % %     change: RF flip angle
% gamma value make it back to 10^6
gammaValue = str2double(get(handles.gammaval, 'String')) * (10^8);
% make it from ms to sec
tauValue   = str2double(get(handles.RFtimeval, 'String')) * (10^-3);
% make it from microtesla to tesla
b1Value    = str2double(get(handles.B1val, 'String')) * (10^-6);

% % % % New calculated value
newThetaValue = gammaValue * b1Value * tauValue;
% % % % Change old theta value with new theta value
% % % %     change it from radians to degrees
% % % %     and take the integer value
set(handles.RFflipangleval, 'String', ...
    num2str(floor(newThetaValue * 180 / pi)));


% --- Executes on button press in updateRFtime.
function updateRFtime_Callback(hObject, eventdata, handles)
% hObject    handle to updateRFtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % % % When RF pulse time (tau) is changed
% % % %     keep constant: B1 magnitude
% % % %     change: RF flip angle
% gamma value make it back to 10^6
gammaValue = str2double(get(handles.gammaval, 'String')) * (10^8);
% make it from ms to sec
tauValue   = str2double(get(handles.RFtimeval, 'String')) * (10^-3);
% make it from microtesla to tesla
b1Value    = str2double(get(handles.B1val, 'String')) * (10^-6);

% % % % New calculated value
newThetaValue = gammaValue * b1Value * tauValue;
% % % % Change old rf angle theta value with new theta value
% % % %     change it from radians to degrees
% % % %     and take the integer value
set(handles.RFflipangleval, 'String', ...
    num2str(floor(newThetaValue * 180 / pi)));


% --- Executes on button press in updateRFangle.
function updateRFangle_Callback(hObject, eventdata, handles)
% hObject    handle to updateRFangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % % % When RF flip angle (theta) is changed
% % % %     keep constant: B1 magnitude
% % % %     change: RF pulse time
% gamma value make it back to 10^6
gammaValue = str2double(get(handles.gammaval, 'String')) * (10^8);
% make it from degrees to radians
thetaValue = str2double(get(handles.RFflipangleval, 'String')) * pi / 180;
% make it from microtesla to tesla
b1Value    = str2double(get(handles.B1val, 'String')) * (10^-6);

% % % % New calculated value
newTauValue = thetaValue / (gammaValue * b1Value);
% % % % Change old tau value with new tau value
% % % %     change it from seconds to milliseconds 
set(handles.RFtimeval, 'String', num2str(newTauValue * (10^3)));


% --- Executes on button press in runSim.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to runSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2



function speedOfSim_Callback(hObject, eventdata, handles)
% hObject    handle to speedOfSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of speedOfSim as text
%        str2double(get(hObject,'String')) returns contents of speedOfSim as a double


% --- Executes during object creation, after setting all properties.
function speedOfSim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speedOfSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
