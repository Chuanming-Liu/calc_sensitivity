% Driver to Calculate Frechet Kernels in terms of Phase Velocity
% NJA, 2014
%
% This involves calling fortran programs plot_wk, frechet, frechet_gv,
% frechet_pv
%
% 2/3/2015 added switch to turn on calculation of psi kernels

function [FRECH_T,FRECH_S] = fmk_kernels_linear(CARD,RUNPATH,PERIODS,ispsi)

% setup_parameters;
% CARD = param.CARDID;
% RUNPATH = param.RUNPATH;
% PERIODS = param.periods;


yaxis = [6100000 6371000];

isfigure = 0;

% Set path to executables
set_plotwk_path;

% Change environment variables to deal with gfortran
setenv('GFORTRAN_STDIN_UNIT', '5')
setenv('GFORTRAN_STDOUT_UNIT', '6')
setenv('GFORTRAN_STDERR_UNIT', '0')

%% Make branch files
disp('--- Make Branch Files ---');
TYPE = 'T';
write_plotwk(TYPE,CARD);

com = ['cat ',RUNPATH,'run_plotwk.t | plot_wk'];
[status,log] = system(com);

TYPE = 'S';
write_plotwk(TYPE,CARD);

com = ['cat ',RUNPATH,'run_plotwk.s | plot_wk'];
[status,log] = system(com);

%% Make Frechet Kernels!
disp('--- Make Frechet Kernels ---');
% set_mineos_path;

TYPE = 'T';
NDISC = 0;
ZDISC = [];
write_frechet(TYPE,CARD,NDISC,ZDISC)

if ispsi == 1
    com = ['cat ',RUNPATH,'run_frechet.t | frechet_psi'];
else
com = ['cat ',RUNPATH,'run_frechet.t | frechet'];
end
[status,log] = system(com);

TYPE = 'S';
write_frechet(TYPE,CARD,NDISC,ZDISC)
% disp('Be patient! This will take ~25 s');
tic
com = ['cat ',RUNPATH,'run_frechet.s | frechet'];
[status,log] = system(com);
toc


%% Make CV Frechet Kernels
disp('--- Make CV Frechet Kernels ---');
branch = 0;

TYPE = 'T';
write_frechcv(TYPE,CARD,branch);

com = ['cat ',RUNPATH,'run_frechcv.t | frechet_cv'];
[status,log] = system(com);

TYPE = 'S';
write_frechcv(TYPE,CARD,branch)

com = ['cat ',RUNPATH,'run_frechcv.s | frechet_cv'];
[status,log] = system(com);

% Convert CV Frechet kernels to ascii
% Will do this for all periods of interest
% Set inside the setparam_MINE.m
disp('--- Convert Frechet CV to ascii ---');

% Program writes run file for draw_frechet_gv, runs it, and reads in
% sensitivity kernels for all periods of interest

TYPE = 'T';
if ispsi == 1
    % psi is in place of sh and sv remains unchanged
    
    FRECH_T = frechpsi_asc(TYPE,CARD,branch,PERIODS);
else
    FRECH_T = frechcv_asc(TYPE,CARD,branch,PERIODS);
end

TYPE = 'S';
FRECH_S = frechcv_asc(TYPE,CARD,branch,PERIODS);

% size(FRECH_T(1).rad)
% size(FRECH_T(1).vsv)
% size(FRECH_T(1).vsh)


if isfigure
    CC = winter(length(PERIODS));
    figure(61)
    clf
    figure(62)
    clf
    xaxis = [0 3E-8];
    for ip = 1:length(PERIODS)
        figure(61)
        subplot(1,2,1)
        hold on
        plot(FRECH_S(ip).vsv,FRECH_S(ip).rad,'-k','linewidth',2,'color',CC(ip,:))
        ylim(yaxis)
        %         xlim(xaxis)
        title('Spheroidal SV')
        subplot(1,2,2)
        hold on
        plot(FRECH_S(ip).vsh,FRECH_S(ip).rad,'--k','linewidth',2,'color',CC(ip,:))
        ylim(yaxis)
        %         xlim(xaxis)
        title('Spheroidal SH')
        
        figure(62)
        subplot(1,2,2)
        hold on
        plot(FRECH_T(ip).vsv,FRECH_T(ip).rad,'-k','linewidth',2,'color',CC(ip,:))
        ylim(yaxis)
        %         xlim(xaxis)
        title('Toroidal SV')
        subplot(1,2,1)
        hold on
        plot(FRECH_T(ip).vsh,FRECH_T(ip).rad,'--k','linewidth',2,'color',CC(ip,:))
        if ispsi == 1
            title('Toroidal PSI')
        else
            title('Toroidal SH')
        end
        ylim(yaxis)
       
    end
%     pause
    %
    %     subplot(1,2,1)
    %     hold on
    %     legend(llegend)
    %     subplot(1,2,2)
    %     hold on
    %     legend(llegend)
end

% savefile = [CARD,'_fcv.mat'];
% save(savefile,'FRECH_T','FRECH_S');
% Change the environment variables back to the way they were
setenv('GFORTRAN_STDIN_UNIT', '-1')
setenv('GFORTRAN_STDOUT_UNIT', '-1')
setenv('GFORTRAN_STDERR_UNIT', '-1')
