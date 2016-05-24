% Write frechet_cv driver
% NJA, 2014
% 
% Must designate mode branch of interest (0 = fundamental)
% 
function write_frechcv(TYPE,CARD,BR)
% setparam_MINE;

% TYPE = 'T';
% CARD = param.CARDID;
% BR = 0; 

setup_parameters;

CARDPATH = param.CARDPATH;
DATAPATH = param.DATAPATH;
TABLEPATH = param.TABLEPATH;
RUNPATH = param.RUNPATH;
SID = param.SID;
TID = param.TID;

if strcmp(TYPE,'T') == 1
%     disp('Toroidal!');
    
    RUNFILE = 'run_frechcv.t';
    TYPEID = TID;
    
elseif strcmp(TYPE,'S') == 1
%     disp('Spheroidal!');
    
    RUNFILE = 'run_frechcv.s';
    TYPEID = SID;
    
else
    disp('No TYPE recognized!');
    
end

if BR == 0
    BRID = '0st';
elseif BR == 1
    BRID = '1st';
elseif BR == 2
    BRID = '2nd';
elseif BR == 3
    BRID = '3rd';
else
    disp('Branch has no name! Change it in the script')
    return
end

% QMOD = [CARDPATH,CARD,'.qmod'];
% BRANCH = [DATAPATH,CARD,'.',TYPEID,'.table_hdr.branch'];
FRECH = [DATAPATH,CARD,'.',TYPEID,'.frech'];
FRECHCV = [DATAPATH,CARD,'.',TYPEID,'.fcv.',BRID];

if exist(FRECHCV,'file') == 2
    % disp('File exists! Removing it now')
    com = ['rm ',FRECHCV];
    [status,log] = system(com);
end

fid = fopen(RUNFILE,'w');
fprintf(fid,'%s\n',FRECH);
fprintf(fid,'%i\n',BR);
fprintf(fid,'%s\n',FRECHCV);
fclose(fid);


com = ['mv run_frechcv* ',RUNPATH,'.'];
[status,log] = system(com);
