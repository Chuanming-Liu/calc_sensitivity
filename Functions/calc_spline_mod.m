% Create spline velocity model from previously defined basis functions and
% spline coefficients
%
% 10/1/2014 NJA
%
% Based off of the last half of spline_S_mc
%
% It's confusing but the majority of hte script has variables named sv, sh,
% ... however, the math is identical to the p-wave case.  To calculate
% spline models for P wave velocities just define the type, pass it p-wave
% spline coefficients and make sure to name the returned values (sv, sh,
% ...) something logical (pv, ph, ...)
%
% 10/7/2014 modified to allow optional argument to turn on capatability
% with psi -- i.e. instead of looking for SV and SH it looks for SV and
% psi.  SV and psi models are output. To creat the SH model just use math!
%
% The variable argument is any vector that isnt empty.  
% Dont get confused because (similar to the treatment of P-wave
% coefficients) the psi values will be stored within vectors named coef_h 
%
% NOTE!! Now the coefficients must be named by coef*.v, coef*.h, coef*.psi
% no more of this coef*.sv, coef*.sh business
%

function [sv,sh,sz,coef] = calc_spline_mod(spline,coefi,TYPE,varargin)


isfigure = 0;

if nargin == 3
    
    % Check to see what wave type we're useing
    if strcmp(TYPE,'S') == 1
%         disp('Generating SV and SH Models');
        coef_v = coefi.v;
        coef_h = coefi.h;
        
    elseif strcmp(TYPE,'P')==1
%         disp('Generating PV and PH Models');
        coef_v = coefi.v;
        coef_h = coefi.h;
    else
        disp('No type found!!')
        return
    end
    
    %% Force Isotropy in the Mantle
    
    % Check to see if models are anisotropic at the bottom
%     if abs(coef_v(end)-coef_h(end))/coef_v(end) > .005
%         %     disp('Forcing Isotropy at the Bottom!')
%         iso_vs = (coef_v(end)+coef_h(end))/2;
%         coef_v(end) = iso_vs;
%         coef_h(end) = iso_vs;
%     end
elseif nargin == 4
%     disp('Generating SV and PSI Models');
    % DONT GET CONFUSED!! FROM HERE ON OUT SH MEANS PSI ...
    coef_v = coefi.v;
    coef_h = coefi.psi;
    
    % Force walk towards isotropy at the bottom of the model
%     iso_vs = (1+coef_h(end))/2;
%     coef_h(end) = iso_vs;
    
    
end

coef_z = coefi.z;


%% Initialize vectors for the splines
z_bs = [spline(1).z,spline(2).z,spline(3).z];
sv_bs = zeros(length(z_bs),1);
sh_bs = zeros(length(z_bs),1);

%% Calculate velocities from splines

iz = 1; % keep count of depth
ic = 1; % keep count of coefficients
for ilayer=1:length(spline)
    clear depth sv_coef sh_coef N
    % Layer specific variables
    [N dum] = size(spline(ilayer).b);
    z = spline(ilayer).z;
    sv_coef = coef_v(ic:N+(ic-1));
    sh_coef = coef_h(ic:N+(ic-1));
    
    for int = 1:length(z)
        
        b = spline(ilayer).b;
        
        for i = 1:N
            sv_bs(iz) = sv_bs(iz)+sv_coef(i)*b(i,int);
            sh_bs(iz) = sh_bs(iz)+sh_coef(i)*b(i,int);
            
        end % end of coefficient loop
        iz = iz+1;
    end % end of local depth loop
    ic = ic+N;
end % end of sed, crust, mant loop

if isfigure
    
    figure(68)
    clf
    hold on
    plot(sv_bs,z_bs,'-','color','r','linewidth',2)
%     plot(sh_bs,z_bs,':','color','r','linewidth',2);
    plot(coef_v,coef_z,'og','linewidth',2)
%     plot(coef_h,coef_z,'ob','linewidth',2);
    set(gca,'ydir','reverse','fontsize',16)
    title(sprintf('TYPE : %s',TYPE));

    figure(69)
    clf
    hold on
    plot(sh_bs,z_bs,':','color','r','linewidth',2);
        
            plot(coef_h,coef_z,'ob','linewidth',2);
            
    set(gca,'ydir','reverse','fontsize',16)
end

% Return SV, SH, and Z
clear sv sh z
sv = sv_bs;
sh = sh_bs;
sz = z_bs;

coef.v = coef_v;
coef.h = coef_h;
coef.z = coef_z;
