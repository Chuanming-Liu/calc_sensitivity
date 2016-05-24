% Check the output from mineos_nohang to make sure no eigenfrequencies got
% stuck.
%
% NJA 8/2014

function check_eigen(TYPE)

setup_parameters;
TABLEPATH = param.TABLEPATH;
MODEPATH = param.MODEPATH;
CARDID = param.CARDID;
TID = param.TID;
% TYPE = 'S';

CARDPATH = pwd;
CARDPATH = [CARDPATH,'/'];

if strcmp(TYPE,'S') == 1
    
    TYPEID = 's0to50';
    
    % Name output files
    ASC = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.asc'];
    
    com = ['tail -n 1 ',ASC];
    [status,result] = system(com);
    
elseif strcmp(TYPE,'T') == 1
    TYPEID = TID;
    RUNFILE = 'run_nohang.t';
    MODENUM = 't.mode';
    
    % Write nohang driver
    ASC = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.asc'];
    
    com = ['tail -n 1 ',ASC];
    [status,result] = system(com);
    
end

if strcmp(result(1:5),'stuck') == 1
    disp(['mineos got stuck at ',FREQ,'!']);
    return
end