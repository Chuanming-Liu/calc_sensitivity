% frechcv_asc.
% Program converts binary frechet kernels in phase velocity to ascii
% versions and saves the output in a mat file
% Involves writing run files for draw_frechet_gv and running the fortran
% program
%
% NJA, 2014

function [FRECH] = frechcv_asc(TYPE,CARD,BRANCH,periods)

%%% Turn on if running not as a function
% TYPE = 'T';
% CARD = param.CARDID;
% BRANCH = 0;

% Get useful info from parameter file
setup_parameters;

CARDPATH = param.CARDPATH;
DATAPATH = param.DATAPATH;
TABLEPATH = param.TABLEPATH;
RUNPATH = param.RUNPATH;
TID = param.TID;
SID = param.SID;
% periods = param.periods;

isfigure = 0;

if strcmp(TYPE,'T') == 1
%     disp('Toroidal!');
    
    RUNFILE = 'run_frechcv_asc.t';
    TYPEID = TID;
    
elseif strcmp(TYPE,'S') == 1
  %   disp('Spheroidal!');
    
    RUNFILE = 'run_frechcv_asc.s';
    TYPEID = SID;
    
else
    disp('No TYPE recognized!');
    
end

if BRANCH == 0
    BRID = '0st';
elseif BRANCH == 1
    BRID = '1st';
elseif BRANCH == 2
    BRID = '2nd';
elseif BRANCH == 3
    BRID = '3rd';
else
    disp('Branch has no name! Change it in the script')
end

FRECHCV = [DATAPATH,CARD,'.',TYPEID,'.fcv.',BRID];




for ip = 1:length(periods)
    
    FRECHASC = [DATAPATH,CARD,'.',TYPEID,'.',BRID,'.',num2str(periods(ip))];
    
    % Write runfile for draw_frechet_gv
    fid = fopen(RUNFILE,'w');
    fprintf(fid,'%s\n',FRECHCV);
    fprintf(fid,'%s\n',FRECHASC);
    fprintf(fid,'%i\n',periods(ip));
    fclose(fid);
    
    % Run draw_frechet_gv
    % disp(sprintf('--- Period : %s',num2str(periods(ip))));
    
    if exist(FRECHASC,'file') == 2
    % disp('File exists! Removing it now')
    com = ['rm ',FRECHASC];
    [status,log] = system(com);
    end
    com = ['cat ',RUNFILE,' | draw_frechet_gv'];
    [status,log] = system(com);
    
    % Load in frechet files for each period
    
    % spheroidal, no aniso: 1=Vs,2=Vp,3=rho
    % spheroidal, aniso: 1=Vsv,2=Vpv,3=Vsh,4=Vph,5=eta,6=rho
    % toroidal, no aniso: 1=Vs,2=rho
    % toroidal, aniso: 1=Vsv,2=Vsh,3=rho
    % disp(FRECHASC)
    fid = fopen(FRECHASC,'r');
    C = textscan(fid,'%f%f%f%f%f%f%f');
    
    if strcmp(TYPE,'S') == 1
        FRECH(ip).per = periods(ip);
        FRECH(ip).rad = C{1};
        FRECH(ip).vsv = C{2};
        FRECH(ip).vpv = C{3};
        FRECH(ip).vsh = C{4};
        FRECH(ip).vph = C{5};
        FRECH(ip).eta = C{6};
        FRECH(ip).rho = C{7};
    elseif strcmp(TYPE,'T') == 1
        FRECH(ip).per = periods(ip);
        FRECH(ip).rad = C{1};
        FRECH(ip).vsv = C{2};
        FRECH(ip).vsh = C{3};
        FRECH(ip).rho = C{4};
    end
    
    if isfigure
        figure(18)
        clf
        plot(FRECH(ip).vsv,(FRECH(ip).rad/1000),'-k','linewidth',2)
        title(sprintf('Per %s',num2str(FRECH(ip).per)));
        ylim([6171 6371])
        % pause
    end
end

% move run files to run directory
com = ['mv run_frechcv_asc.? ',RUNPATH,'.'];
[status,log] = system(com);
