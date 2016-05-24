% Setup Parameters for Abridged Version of Monte Carlo Inversion for Radial
% Anisotropy
%
% For running on Natalies mac mini
%
% NJA, 9/30/2014
%
% Based off shear inversion codes from Ge Jin
% psi controls the ratio of SV to SH


% ----------- Parameters for Array ----------- 
param.lalim=[6 14.5];
param.lolim=[36 44];
param.gridsize=0.5;   % in degrees (larger setting for broader array)
param.component='LH';

%  ----------- Parameters for Splines ----------- 
param.sed_spline = [1];
param.crust_spline = [2 3 4 5]; % index of splines within the crust

% ----------- Parameters from Phase Velocities ----------- 
% param.Cperiods = [25 32 40 50 60 80 100];
% param.LAperiods = [16.667 13.3333 11.1111 9.5238 8.3333];
% param.RAperiods = [22.2222 16.6667 13.3333 11.1111 9.5238 8.3333];
% param.allperiods = [param.Cperiods,param.RAperiods];
% param.Uperiods = [22.222 16.6667 13.333 11.1111 9.5238 8.333];

% ----------- Parameters for Grid Search ----------- 
% param.bot = [210 230];
param.bot = 300;
param.sedh = [1 3];
param.rho = [2570 2967 3330];
param.sedv = [2500];
param.crusth_ratio = [1.0];
param.crustv = [3400 3500];
param.crustpsi = [1];
param.mantlev = [3900 4000 4100];
param.mantlepsi = [0.98 1 1.05];
param.mantle_structS = [1.0266,1.0082,1,1.0123,1.0294;1,1,1,1,1;0.974,0.9919,1,0.988,0.986];
param.mantle_structP = [0.9917,0.9808,1,1.0072,1.0209;1,1,1,1,1;1.007,1.0150,1,0.993,0.982];
param.psi_struct = [1,1,1,1,1]; % scalings for bulk Vsh/Vsv
param.vpvs = 1.8;
param.maxchi2 = 2.5;

% ----------- Parameters for Running Mineos ----------- 
param.CARDID = 'goc_nolake5';
% param.CARDID = 'goc_lake4';
param.TID = 't0to150';
param.SID = 's0to150';
% param.TABLEPATH = '/home/segment/Natalie/MODES/Data/'; % MCS
% param.MODEPATH = '/home/segment/Natalie/MODES/mk_modes/MODE.in/'; %MCS
param.TABLEPATH = '/Users/naccardo/Unix/MINEOS/MODE_tables/';
param.INPUTPATH = '/Users/naccardo/Dropbox/Unix_DB/mcmc/data/';
param.MODEPATH = param.INPUTPATH;
% param.MODEPATH = '/Users/naccardo/Unix/MINEOS/MODE_tables/MODE.in/';
param.DATAPATH = [param.TABLEPATH,param.CARDID,'/tables/'];
if ~exist(param.DATAPATH,'dir')
    mkdir(param.DATAPATH)
end
param.RUNSPATH = [param.TABLEPATH,param.CARDID,'/runs/'];
param.CARDPATH = pwd;
param.CARDPATH = [param.CARDPATH,'/'];
param.RUNPATH = pwd;
param.RUNPATH = [param.RUNPATH,'/run/'];

% ----------- parameters for inversion ----------- 
param.maxiter = 4;
param.maxpert = 100; % amount we will let spline coef change (m/s)

% ----------- Parameters for Mineos Programs ----------- 
% param.eps = 1e-15;
% param.wgrav = 1000;
% param.jcom = 2;
% param.lmin = 0;
% param.lmax = 20000;
% param.nmin = 0;
% param.nmax = 0;
param.maxN = 18000; % max number of modes
param.minF = 0; % min frequency in mHz -- should match file names
param.maxF = 150.05; % max frequency in mHz
param.minL = 0; % min angular order
param.maxL = 6000; % max angular order
param.wmin = 0.05; % should match minF and maxF
param.wmax = 151;


% ----------- Parameters for Monte Carlo Inversion ----------- 
param.testN = 60;
param.crusth_var = 5;
param.sedh_var = 1;
param.velocity_var = 0.065;
param.psi_var = 0.035;
param.r = 0.2;

% parameters for plotting
LBLFONT = 17;
LBLFONTNAME = 'Times New Roman';

param.data = '/Users/naccardo/Dropbox/Unix_DB/mcmc/progs/surf96_mcmc/data/LR_phasev_results.mat';
