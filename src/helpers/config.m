% % % % IRINA GRIGORESCU
% % % % DATE: 13-Feb-2017
% % % % Config parameters
% % % % 

function params = config()
% % % Just runs some stuff needed to be stored in memory 
% % % of a different program

%% The parameters which will be passed as output
params = struct;

%% Tissue structure with spin density, T1 and T2 
params.tissue = struct;

params.tissue.csf.rho =    1.0;
params.tissue.csf.t1  = 4500;
params.tissue.csf.t2  = 2200;

params.tissue.wm.rho =   0.65;
params.tissue.wm.t1  = 600;
params.tissue.wm.t2  =  80;

params.tissue.gm.rho =   0.8;
params.tissue.gm.t1  = 950;
params.tissue.gm.t2  = 100;

params.tissue.fat.rho =  0.9;
params.tissue.fat.t1 = 250;
params.tissue.fat.t2 =  60;


end
