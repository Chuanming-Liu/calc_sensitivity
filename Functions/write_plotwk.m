% Write driver for plot_wk
% NJA, 2014

function write_plotwk(TYPE,CARD)

setup_parameters;

DATAPATH = param.DATAPATH;
RUNPATH = param.RUNPATH;
TID = param.TID;
SID = param.SID;

% Parameters for searching for modes
% Wmin = 0;
% Wmax = 51;
wmin = param.wmin;
wmax = param.wmax;

% TYPE = 'T';
% CARD = param.CARDID;

if strcmp(TYPE,'T') == 1

    RUNFILE = 'run_plotwk.t';
    
    TABLEHDR = [DATAPATH,CARD,'.',TID,'.table_hdr'];
    

elseif strcmp(TYPE,'S') == 1
    
    RUNFILE = 'run_plotwk.s';
    TABLEHDR = [DATAPATH,CARD,'.',SID,'.table_hdr'];
    
else
    
    disp('Type does not exist!')

end

% Write out instructions to run file
fid = fopen(RUNFILE,'w');
fprintf(fid,'table %s\n',TABLEHDR);
fprintf(fid,'search \n');
fprintf(fid,'1 %4.2f %4.2f \n',[wmin wmax]);
fprintf(fid,'99 0 0 \n');
fprintf(fid,'branch \n');
fprintf(fid,'\n');
fprintf(fid,'quit \n');
fclose(fid);

% Move to run directory

com = ['mv run_plotwk* ',RUNPATH,'.'];
[status,log] = system(com);

