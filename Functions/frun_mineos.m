% Function version of Run Minos Bran
% Original version is name run_minos.m
%
%
% NJA, 2014
%
% Changed 8/2014 to use only one run file for spheroidal (instead of 6
% separate)
%
% 9/2014 modified to run codes out of the run directory
%
% 10/27/2014 Modified to return status of nohang commands
%
% 1/30/2015 modifed to use q corrected phae velocities and allow choice to
% run mineos strip and table

function [status_S,status_T] = frun_mineos(CARD,runtable)

setup_parameters;
RUNPATH = param.RUNPATH;

% % Turn on if only want to calculate S or T
% SONLY = 0;
% TONLY = 1;

% Check to see if run directory exists and create it if not

if exist('run','dir') ~= 7
    mkdir('run');
end

% Set path to executables
% set_mineos_path;

% Change environment variables to deal with gfortran
setenv('GFORTRAN_STDIN_UNIT', '5') 
setenv('GFORTRAN_STDOUT_UNIT', '6') 
setenv('GFORTRAN_STDERR_UNIT', '0')


%% Run Spheroidal Branches First
disp('---- Calculate Spheroidal Mode Set 1 ----')
TYPE = 'S';
% Write out run files -- be sure that paths are correct!
write_mineos_drivers(TYPE,CARD);

% mineos_nohang for s1-s5
% disp('Running mineos_nohang');
com = ['cat ',RUNPATH,'run_nohang.s | mineos_nohang'];
% tic
[status_S(1),log] = system(com);


% disp('---- Calculate Spheroidal Mode Set 2 ----')
% com = ['cat ',RUNPATH,'run_nohang2.s | mineos_nohang'];
% [status_S(2),log] = system(com);
% toc

% Check to make sure nothing got hung up
% disp('Checking for Stuck Eigenfrequencies')

check_eigen('S');

% % mineos_q
% disp('Running mineos_q');
% com = ['cat ',RUNPATH,'run_q.s | mineos_q'];
com = ['cat ',RUNPATH,'run_q.s | mineos_qcorrectphv'];
[status,log] = system(com);
% 
if runtable == 1
% mineos_strip
disp('Running mineos_strip');
com = ['cat ',RUNPATH,'run_strip.s | mineos_strip'];
[status,log] = system(com);
% 
% % mineos_table
% disp('Running mineos_table');
com = ['cat ',RUNPATH,'run_table.s | mineos_table'];
[status,log] = system(com);
end

%% Run toroidal branches next
disp('---- Calculate Toroidal Modes ----')
TYPE = 'T';
% Write out run files -- be sure that paths are correct!
write_mineos_drivers(TYPE,CARD);

% mineos_nohang
% disp('Running mineos_nohang');

com = ['cat ',RUNPATH,'run_nohang.t | mineos_nohang'];
[status_T,log] = system(com);

% Check to make sure nothing got hung up
% disp('Checking for Stuck Eigenfrequencies')

check_eigen('T');

% % mineos_q
% disp('Running mineos_q');
% com = ['cat ',RUNPATH,'run_q.t | mineos_q'];
com = ['cat ',RUNPATH,'run_q.t | mineos_qcorrectphv'];
[status,log] = system(com);
% 
if runtable == 1
% mineos_strip
disp('Running mineos_strip');
com = ['cat ',RUNPATH,'run_strip.t | mineos_strip'];
[status,log] = system(com);
% 
% % mineos_table
% disp('Running mineos_table');
com = ['cat ',RUNPATH,'run_table.t | mineos_table'];
[status,log] = system(com);
end
% Change the environment variables back to the way they were
setenv('GFORTRAN_STDIN_UNIT', '-1') 
setenv('GFORTRAN_STDOUT_UNIT', '-1') 
setenv('GFORTRAN_STDERR_UNIT', '-1')
