% Inversion code to invert phase velocities for perturbations to spline
% coefficients
%
% NJA, 11/22/2014
%
% Based off of code from 5/2014
%

function mest = finverse_aniso(Lphv,Lstd,Rphv,Rstd,forward,kern)


isfigure = 1;

setup_parameters;
periods = param.Cperiods;
% Get forward model
tphv = forward.tphv;
tper = forward.tper;
sphv = forward.sphv;
sper = forward.sper;

for ip = 1:length(periods)
    % Double check that periods match
    
    tmatch = abs(tper - periods(ip));
    tind = find(tmatch == min(tmatch));
    smatch = abs(sper - periods(ip));
    sind = find(smatch == min(smatch));
	% disp(tper)
	% length(tphv(tind))
	% length(sphv(sind))
    
    dCT(ip) = (Lphv(ip)-tphv(tind));
    dCS(ip) = (Rphv(ip)-sphv(sind));
end

dobs = [dCT,dCS]';

% change to m/s
dobs = dobs*1000;

disp(dobs)

%% Build G matrix - pull in coefficient specific sensitivities from the
% perturb_mod script
% This will be organized with the Love wave sensitivities first and then
% the VSV sensitivities

dTSVdCOEF = kern.TSV;
dTSHdCOEF = kern.TSH;
dSSVdCOEF = kern.SSV;
dSSHdCOEF = kern.SSH;

% check to make sure matrices look correct
if size(dTSVdCOEF) ~= size(dSSVdCOEF)
    disp('Error! Kernel Matrices dont match!')
    error;
elseif size(dTSHdCOEF) ~= size(dSSHdCOEF)
        disp('Error! Kernel Matrices dont match!')
    error;
end

[Ncoef Nper] = size(dTSHdCOEF);



% Firt SV
GSV = zeros(Nper*2,Ncoef);
for icoef = 1:Ncoef
    for ip = 1:Nper
        % Toroidal first
        GSV(ip,icoef) = dTSVdCOEF(icoef,ip);
        % Then spheroidal
        GSV(ip+Nper,icoef) = dSSVdCOEF(icoef,ip);
    end
end

% Next SH

GSH = zeros(Nper*2,Ncoef);
for icoef = 1:Ncoef
    for ip = 1:Nper
        % Toroidal first
        GSH(ip,icoef) = dTSHdCOEF(icoef,ip);
        % Then spheroidal
%         GSH(ip+Nper,icoef) = dSSHdCOEF(icoef,ip);
GSH(ip+Nper,icoef) = 0;
    end
end

if isfigure
    % Lets check the G matrix
figure(1)
clf
CC=jet(Ncoef);
for icoef = 1:Ncoef
    hold on
    plot(periods,GSV(1:7,icoef),'-','color',CC(icoef,:),'linewidth',2)
    plot(periods,GSV(8:end,icoef),':','color',CC(icoef,:),'linewidth',2)
% pause
end
ylabel('dAdC')
xlabel('Periods (s)')
title('Sensitivities to SV')
hleg = legend('Toroidal','Spheroidal','location','northeast');

figure(2)
clf
CC=jet(Ncoef);
for icoef = 1:Ncoef
    hold on
    plot(periods,GSH(1:7,icoef),'-','color',CC(icoef,:),'linewidth',2)
    plot(periods,GSH(8:end,icoef),':','color',CC(icoef,:),'linewidth',2)
% pause
end
ylabel('dAdC')
xlabel('Periods (s)')
title('Sensitivities to SH')
hleg = legend('Toroidal','Spheroidal','location','northeast');


% Lets compare the observed and estimated phase velocities
figure(3)
clf
hold on
plot(periods,Lphv,'-k','linewidth',2)
plot(periods,Rphv,'-','color','r','linewidth',2)
plot(periods,tphv,'ok','linewidth',2)
plot(periods,sphv,'or','linewidth',2)
xlabel('Periods')
ylabel('Phase Velocity')
leg = legend('Obs L Phv','Obs R Phv','Est T Phv','Est S Phv');
end
% 
% figure(4)
% clf
% hold on
% plot([periods,periods],sum(GSV,2),'or','linewidth',2)
% plot([periods,periods],sum(GSH,2),'ok','linewidth',2)

%% Genearte the weight matrix -- based off of ratios of standard deviation
% Most accurate measurements will have 

% Inverse the matrix to get mest (perturbations to the coeffi
% mest.SV = (GSV'*GSV)\(GSV'*dobs);

mest.SV = lsqr(GSV,dobs);

GSV*mest.SV
GSV

% mest.SV = mest.SV/(1000^3);
mest.SV = mest.SV/(10^8);

disp('SV')
disp(mest.SV);

% mest.SH = (GSH'*GSH)\(GSH'*dobs);
mest.SH = lsqr(GSH,dobs);

mest.SH = mest.SH/(10^8);
% mest.SH = mest.SH/(1000^3);

% disp('SH')
% disp(mest.SH);
% 
% save(savefile,'mest');

%     end % Ny
% end % Nx
