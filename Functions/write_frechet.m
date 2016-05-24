% Write run file for frechet.f
% NJA, 2014
%
% 8/2014 Modified to run with only one spheroidal run flie
function write_frechet(TYPE,CARD,NDISC,ZDISC)
setup_parameters;

% TYPE = 'S';
% CARD = param.CARDID;
% NDISC = 0;
% ZDISC = [];


CARDPATH = param.CARDPATH;
DATAPATH = param.DATAPATH;
TABLEPATH = param.TABLEPATH;
RUNPATH = param.RUNPATH;
TID = param.TID;
SID = param.SID;

% Set number of discontinuities to add
if NDISC > 0
    
    if length(ZDISC) > NDISC
        disp('Mismatch in discontinuity depths!')
    end
end

if strcmp(TYPE,'T') == 1
%    disp('Toroidal!');
    
    RUNFILE = 'run_frechet.t';
    TYPEID = TID;
    EIG = [DATAPATH,CARD,'.',TYPEID,'.eig'];
    
elseif strcmp(TYPE,'S') == 1
 %   disp('Spheroidal!');
    
    RUNFILE = 'run_frechet.s';
    TYPEID = SID;
    EIG = [DATAPATH,CARD,'.',TYPEID,'.eig'];
else
    disp('No TYPE recognized!');
end

QMOD = [CARDPATH,CARD,'.qmod'];
BRANCH = [DATAPATH,CARD,'.',TYPEID,'.table_hdr.branch'];
FRECH = [DATAPATH,CARD,'.',TYPEID,'.frech'];

% Check to see if file exists ... program will not overwrite it if it does
if exist(FRECH,'file') == 2
    % disp('Frechet file already exists! Removing it now')
    com = ['rm ',FRECH];
    [status,log] = system(com);
end

% Write runfile
fid = fopen(RUNFILE,'w');
fprintf(fid,'%s\n',QMOD);
fprintf(fid,'%s\n',BRANCH);
fprintf(fid,'%s\n',FRECH);
fprintf(fid,'%s\n',EIG);

fprintf(fid,'%i\n',NDISC);

for idisc = 1:NDISC
    fprintf(fid,'%i\n',ZDISC(idisc));
end
fprintf(fid,'\n');
fclose(fid);

% Move run files to run directory
com = ['mv run_frechet.? ',RUNPATH,'.'];
[status,log] = system(com);
