% Test the predicted phase velocity perturbations due to our two different
% parameterizations
%
% Two roads diverged in a yellow wood,
% And sorry I could not travel both
%
% NJA, 2/3/2015

clear

isfigure = 0;
savefile = 'test_kerns.mat';
ischeckmod = 0;

%% Setup parameters
setup_parameters;

vpvs = param.vpvs;
bot = param.bot;
sedh = param.sedh;
sedv = param.sedv;
crusth_ratio = param.crusth_ratio;
crustv = param.crustv;
mantlev = param.mantlev;
mantlepsi = param.mantlepsi;
mantstructS = param.mantle_structS;
mantstructP = param.mantle_structP;
mantstructPSI = param.psi_struct;
component = param.component;
Tperiods = param.Cperiods;
RAperiods = param.RAperiods;
LAperiods = param.LAperiods;
runpath = param.RUNPATH;
Allperiods = [Tperiods,RAperiods];
% cardid = param.CARDID;

icardid = 'simple';

% Parameters for inversion
% maxiter = param.maxiter;
maxiter = 10;
maxpert = param.maxpert;
maxpert = 100;

% Parameters to loop through the array
lalim = param.lalim;
lolim = param.lolim;
gridsize = param.gridsize;
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);
[xi yi]=ndgrid(xnode,ynode);

crusth = 5;
sedh = 17;

% % Choose parameters for the starting model
% sedh1 = sedh(2);
% sedv1 = sedv(1);
% crustv1 = crustv(2);
% mantlev1 = mantlev(2);
% mantlepsi1 = mantlepsi(3);
% mantstructS1 = mantstructS(1,:);
% mantstructP1 = mantstructP(1,:);

% parameters from phase velocity
periods = param.Cperiods;

% turn on use of psi parameterization
ispsi = 1;
ischeckmodel = 0;

maxchi2 = param.maxchi2;

% Set path to executables
set_mineos_path;
% addpath('/home/segment/Natalie/bspline');

% Load phase velocity information
ddata = param.data;
load(ddata);

%% Start Inversion
ilat = 5;
ilon = 5;

count = 0;

% Get initial model from goc card

% lmodel = read_model_card([icardid,'.card']);
%
% radind = find(lmodel.rad > (6371-bot)*1000);
% ivsv = lmodel.vsv(radind);
% ivsh = lmodel.vsh(radind);
% ivpv = lmodel.vpv(radind);
% ivph = lmodel.vph(radind);
% irad = lmodel.rad(radind);
% cvsv = ivsv;
% cvsh = ivsh;
% cvpv = ivpv;
% cvph = ivph;

% initialize important vectors
iter = 1;
count = count+1;
ilon = ilon+1;

disp(['PROCESSING ',num2str(xi(ilat,ilon)),'x',num2str(yi(ilat,ilon))]);
% Get grid point observations

% ------ Get observations -----
% Teleseismic
for ip=1:length(Tperiods)
    if ~isempty(phasev(ilat,ilon).RTphv)
        obs.RTphv(ip) = phasev(ilat,ilon).RTphv(ip);
        obs.RTstd(ip) = phasev(ilat,ilon).RTstd(ip);
        
    else
        obs.RTphv(ip) = NaN;
        obs.RTstd(ip) = NaN;
    end
    if ~isempty(phasev(ilat,ilon).LTphv)
        obs.LTphv(ip) = phasev(ilat,ilon).LTphv(ip);
        obs.LTstd(ip) = phasev(ilat,ilon).LTstd(ip);
    else
        obs.LTphv(ip) = NaN;
        obs.LTstd(ip) = NaN;
        
    end
    
end

% Find locations with weird std
Rzero = find(obs.RTstd == 0);
Lzero = find(obs.LTstd == 0);

obs.LTstd(Lzero) = 0.1;
obs.RTstd(Rzero) = 0.1;
nanR = find(isnan(obs.RTphv)==1);
nanL = find(isnan(obs.LTphv)==1);
% Ambient Noise
for ip = 1:length(LAperiods)
    if ~isempty(phasev(ilat,ilon).LAphv)
        obs.LAphv(ip) = phasev(ilat,ilon).LAphv(ip);
        %                 obs.LAerr(ip) = phasev(ilat,ilon).LAerr(ip);
        obs.LAerr(ip) = 0.15;
    else
        obs.LAphv(ip) = NaN;
        obs.LAerr(ip) = NaN;
    end
end
for ip = 1:length(RAperiods)
    if ~isempty(phasev(ilat,ilon).RAphv)
        obs.RAphv(ip) = phasev(ilat,ilon).RAphv(ip);
        %                 obs.RAerr(ip) = phasev(ilat,ilon).RAerr(ip);
        obs.RAerr(ip) = 0.15;
    else
        obs.RAphv(ip) = NaN;
        obs.RAerr(ip) = NaN;
    end
end


periods = Allperiods;
obs.Rphv = [obs.RTphv,obs.RAphv];
obs.Lphv = [obs.LTphv,NaN,obs.LAphv];
obs.Lstd = [obs.LTstd,NaN,obs.LAerr];
obs.Rstd = [obs.RTstd,obs.RAerr];

% Check to make sure things aren't all NaN
if length(nanL) > 2 || length(nanR) > 2
    disp('All NaN data! Going to next point')
    return
end

if exist(savefile) ~= 2
    % Create the initial model --
    
    
    sedz = [0:5];
    crustz = [5:30];
%     mantz = [30:20:150];
    smantz = [30:5:390];
%     lmantz = [131:390];
    
    depz = [sedz,crustz,smantz];
    
    vsv = [2500*ones(size(sedz)),3500*ones(size(crustz)),4200*ones(size(smantz))];
    vsh = vsv;
    vpv = vsv*1.82;
    vph = vpv;
    
    psi = (vsh.^2)./(vsv.^2);
    
    write_MINE_mod(smantz(end),sedz(end),crustz(end),vsv',vsh',depz',vpv',vph',['simple.card']);
    
    if ischeckmod == 1
        check_model('simple.card')
    end
    %% Create kernels
    % first run mineos once
    %         frun_mineos([icardid,'.card']);
    
    %         return
    
    cardid = 'simple';
    frun_mineos([cardid,'.card'],1);
    
    % -------  Create Psi & SV Kernels -------
    [frech_Tpsi,frech_Spsi] = fmk_kernels_linear(icardid,runpath,periods,1);
    pause
    % ------- Create SV and SH kernels -------
    [frech_T,frech_S] = fmk_kernels_linear(icardid,runpath,periods,0);
    
    
    % Integrate over depth to calculate perturbation to phase velocity
    
    [qcor,ucor] = calc_qphv(periods);
    iforward.tphv = qcor.tphv;
    iforward.tper = qcor.tper;
    iforward.sphv = qcor.sphv;
    iforward.sper = qcor.sper;
    
    
    % Perturb the initial model by 1% and run mineos again
    nvsv = vsv+(0.01*vsv);
    nvsh = vsh+(0.02*vsh);
    nvpv = vpv+(0.01*vpv);
    nvph = vph+(0.02*vph);
    
    npsi = (nvsh.^2)./(nvsv.^2);
    
    write_MINE_mod(smantz(end),sedz(end),crustz(end),nvsv',nvsh',depz',nvpv',nvph',['nsimple.card']);
    
    cardid = 'nsimple';
    frun_mineos([cardid,'.card'],0);
    
    clear qcor ucor
    [qcor,ucor] = calc_qphv(periods);
    nforward.tphv = qcor.tphv;
    nforward.tper = qcor.tper;
    nforward.sphv = qcor.sphv;
    nforward.sper = qcor.sper;
    
    save(savefile);
else
    load(savefile)
end
if isfigure
figure(1)
clf
subplot(1,2,1)
hold on
plot(vsv-vsh,depz,'-k','linewidth',1)
plot(nvsv-nvsh,depz,'-g','linewidth',1)
subplot(1,2,2)
hold on
plot(vpv-vph,depz,'-k','linewidth',1)
plot(nvpv-nvph,depz,'-g','linewidth',1)
end
figure(2)
clf
subplot(1,2,1)
hold on
plot(iforward.sper,iforward.sphv,'ob','linewidth',2)
plot(nforward.sper,nforward.sphv,'xb','linewidth',2)
xlim([5 100])
title('Spheroidal Phv')
subplot(1,2,2)
hold on
plot(nforward.tper,nforward.tphv,'xr','linewidth',2)
plot(iforward.tper,iforward.tphv,'or','linewidth',2)
xlim([5 100])
title('Toroidal Phv')

% Calculate actual perturbations to phase velocity
obsdT = nforward.tphv-iforward.tphv;
obsdS = nforward.sphv-iforward.sphv;

figure(3)
clf
subplot(1,2,1)
hold on
% plot(iforward.sper,obsdS,'xk','linewidth',2)
title('dC Spheroidal')
xlim([5 100])
subplot(1,2,2)
hold on
% plot(iforward.tper,obsdT,'xk','linewidth',2)
title('dC Toroidal')
xlim([5 100])

if isfigure
CC = spring(length(periods));
DC = winter(length(periods));
yaxis = [6100000 6371000];
figure(61)
clf
figure(63)
clf
figure(62)
clf
figure(64)
clf
for ip = 1:length(periods)
figure(61)
subplot(1,2,1)
hold on
plot(frech_S(ip).vsv,frech_S(ip).rad,'-k','linewidth',2,'color',CC(ip,:))
ylim(yaxis)
%         xlim(xaxis)
title('Spheroidal SV')
subplot(1,2,2)
hold on
plot(frech_S(ip).vsh,frech_S(ip).rad,'--k','linewidth',2,'color',CC(ip,:))
ylim(yaxis)
%         xlim(xaxis)
title('Spheroidal SH')

figure(62)
subplot(1,2,2)
hold on
plot(frech_T(ip).vsv,frech_T(ip).rad,'-k','linewidth',2,'color',CC(ip,:))
ylim(yaxis)
%         xlim(xaxis)
title('Toroidal SV')
subplot(1,2,1)
hold on
plot(frech_T(ip).vsh,frech_T(ip).rad,'--k','linewidth',2,'color',CC(ip,:))
ylim(yaxis)


figure(63)
subplot(1,2,1)
hold on
plot(frech_Spsi(ip).vsv,frech_Spsi(ip).rad,'-k','linewidth',2,'color',DC(ip,:))
ylim(yaxis)
%         xlim(xaxis)
title('Spheroidal SV')
subplot(1,2,2)
hold on
plot(frech_Spsi(ip).vsh,frech_Spsi(ip).rad,'--k','linewidth',2,'color',DC(ip,:))
ylim(yaxis)
%         xlim(xaxis)
title('Spheroidal SH')

figure(64)
subplot(1,2,2)
hold on
plot(frech_Tpsi(ip).vsv,frech_Tpsi(ip).rad,'-k','linewidth',2,'color',DC(ip,:))
ylim(yaxis)
%         xlim(xaxis)
title('Toroidal SV')
subplot(1,2,1)
hold on
plot(frech_Tpsi(ip).vsh,frech_Tpsi(ip).rad,'--k','linewidth',2,'color',DC(ip,:))
ylim(yaxis)
title('Toroidal PSI')


% disp(ip)

end
end
% Calculate perturbations via integration over
N = length(frech_T(1).rad);

for ip = 1:length(periods)
    clear dCT1 dCS1 dCT2 dCS2
    if ip == 1
        figure(7)
        clf
    end
    for iz = 1:length(depz)-1
        dvsv = nvsv(iz)-vsv(iz);
        dvsh = nvsh(iz)-vsh(iz);
        dpsi = npsi(iz)-psi(iz);
        dz = abs(depz(iz)-depz(iz+1))*1000;
        trad = length(frech_T(ip).rad)-(iz-1);
        srad = length(frech_S(ip).rad)-(iz-1);
        
        if depz(iz)-(frech_S(ip).rad(srad))>0
            error('Depths do not match for S!')
        elseif depz(iz)-(frech_T(ip).rad(trad))>0
            error('Depths do not match for T')
        end
%         figure(7)
%         hold on
%         plot(frech_S(ip).vsv(radind),6371-frech_S(ip).rad(radind)/1000,'og','linewidth',2)
        % calculate via sv and sh parameterization
        dCT1(iz) = (frech_T(ip).vsh(trad)*dvsh+frech_T(ip).vsv(trad)*dvsv)*dz;
        dCS1(iz) = (frech_S(ip).vsh(srad)*dvsh+frech_S(ip).vsv(srad)*dvsv)*dz;
        % calculate via psi and sv parameterization
%         dCT2(iz) = ((frech_Tpsi(ip).vsh(trad)*(nvsv(iz)^2)/(2*nvsh(iz))*dpsi)+(frech_Tpsi(ip).vsv(trad)*dvsv))*dz;
        dCT2(iz) = ((frech_Tpsi(ip).vsh(trad)*dpsi)+(frech_Tpsi(ip).vsv(trad)*dvsv))*dz;
        dCS2(iz) = (frech_Spsi(ip).vsv(srad)*dvsv)*dz;
    end
    
    pertCT1(ip) = sum(dCT1);
    pertCS1(ip) = sum(dCS1);
    pertCT2(ip) = sum(dCT2);
    pertCS2(ip) = sum(dCS2);
end


disp('Perturbations to Love')
disp(pertCT1)
disp('Perturbations to Rayleigh')
disp(pertCS1)
% convert to km/s



figure(3)
subplot(1,2,1)
hold on
plot(periods,pertCS1(ip),'ob','linewidth',2)
plot(periods,pertCS2(ip),'og','linewidth',2)
subplot(1,2,2)
hold on
plot(periods,pertCT1,'ob','linewidth',2)
plot(periods,pertCT2,'og','linewidth',2)

figure(4)
clf
subplot(1,2,1)
hold on
plot(periods,(pertCS1(ip)-pertCS2(ip))./(pertCS1)*100,'ob','linewidth',2)
ylabel('Percent Difference')
title('Spheroidal')
subplot(1,2,2)
hold on
plot(periods,(pertCT1(ip)-pertCT2(ip))./pertCT1*100,'or','linewidth',2)
ylabel('Percent Difference')
title('Toroidal')