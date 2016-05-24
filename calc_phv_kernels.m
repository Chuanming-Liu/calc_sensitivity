% Create new velocity models and calculate the resulting sensitivity
% kernels for isotropic models
%
% NJA, 4/20/2016

clear

%% Setup parameters
setup_parameters;
addpath('./Functions/')
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
% periods = param.Uperiods;
runpath = param.RUNPATH;
cardid = param.CARDID;

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

isfigure = 1;

% parameters from phase velocity
% periods = param.Cperiods;

ANperiods = [7.68 8.53 9.58 10.93 12.72 15.21 18.9];
periods = [ANperiods,25,32,40,50,60,80];
% turn on use of psi parameterization
ispsi = 1;
ischeckmodel = 0;

maxchi2 = param.maxchi2;

% Set path to executables
set_mineos_path;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Set up parameters for new model card  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(strfind(cardid,'_nolake'))
    islake = 1;
else
    islake = 0;
end

%%% Remake the velocity model? If yes then mineos will run
makemodel = 0;
if makemodel == 1;
    isrestart = 1;
else
    
isrestart = 0;
end

if makemodel
    
    %%% Only parameters to change!!
    vecdisc = [2 35 100]; % layers - stop at the mantle
    veczkm = [1 4 10];
    vecptop = [1.8 5 7.6]; % P velocity at the shallow boundary
    vecpbot = [1.8 7 8.4]; % P velocity at the deep boundary
    
    
    vecstop = vecptop/1.85;
    vecsbot = vecpbot/1.85;
    vecrho = [1.5 2.9 3.38 ]; % g/m^3
    vecQK = [50000.0  50000.0 50000.0];
    vecQS = [100.0 300.0 300.0];
    
    if islake
        % Add lake layer if we want a lake
        vecdisc = [0.6 vecdisc];
        veczkm = [0.6 veczkm];
        vecptop = [1.5 vecptop];
        vecpbot = [1.5 vecpbot];
        vecstop = vecptop/1.85;
        vecstop(1) = 0;
        vecsbot = vecpbot/1.85;
        vecrho = [1 vecrho];
        vecQK = [57823.0 vecQK];
        vecQS = [999999.0 vecQS];
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Convert model to depth %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nlayer = numel(vecdisc);
for ii = 1:Nlayer-1
    % Calculate the slope for the parameter
    mP = (vecpbot(ii)-vecptop(ii))/vecdisc(ii);
    mS = (vecsbot(ii)-vecstop(ii))/vecdisc(ii);
    mR = 0;
    mQK = 0;
    mQS = 0;
    % Figure out how many slices per layer
    step = veczkm(ii);
    Nstep = floor(vecdisc(ii)/step);
    if ii == 1
        Z = 0;
        P = vecptop(ii);
        S = vecstop(ii);
        R = vecrho(ii);
        QK = vecQK(ii);
        QS = vecQS(ii);
    end
    Zbot = Z(end)+vecdisc(ii);
    count = 0;
    
    while Z(end) <= Zbot
        count = count+1;
        iz = Z(end)+step;
        if iz > Zbot
                P = [P vecptop(ii+1)];
                S = [S vecstop(ii+1)];
                R = [R vecrho(ii+1)];
                QK = [QK vecQK(ii+1)];
                QS = [QS vecQS(ii+1)];
                Z = [Z Z(end)];
            break
        end
        ip = vecptop(ii)+mP*iz;
        is = vecstop(ii)+mS*iz;
        ir = vecrho(ii)+mR*iz;
        iQK = vecQK(ii)+mQK*iz;
        iQS = vecQS(ii)+mQS*iz;
        
        P = [P ip];
        S = [S is];
        R = [R ir];
        QK = [QK iQK];
        QS = [QS iQS];
        Z = [Z iz];
    end
    
end % ii

for ii = 1:length(Z)-1
    tempZ(ii) = Z(end-(ii-1));
    tempS(ii) = S(end-(ii-1));
    tempP(ii) = P(end-(ii-1));
    tempQK(ii) = QK(end-(ii-1));
    tempQS(ii) = QS(end-(ii-1));
    tempR(ii) = R(end-(ii-1));
end
P = [tempP,vecptop(1)];
QK = [tempQK,vecQK(1)];
QS = [tempQS,vecQS(1)];
R = [tempR,vecrho(1)];
if islake
    S = [tempS, 0];
    S(end-1) = 0;
else
    S = [tempS, vecstop(1)];
end
Z = [tempZ,0];



if isfigure
figure(1)
clf
hold on
plot(P,Z,'-o')
plot(S,Z,'-o','color','r')
plot(R,Z,'-o','color','g')
set(gca,'ydir','reverse')
grid on
ylim([0 45])

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Write out model card %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    create_Mineos_cards(S'*1e3,Z',P'*1e3,QK',QS',R'*1e3,[cardid,'.card']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Write out the qmodel %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% match qmodel only based on presence/absence of lake
if vecrho(1) == 1
    system(['cp goc_lake.qmod ',cardid,'.qmod']);
else
    system(['cp goc.qmod ',cardid,'.qmod']);
end
end % makemodel
% Get initial model from goc card

lmodel = read_model_card([cardid,'.card']);

radind = find(lmodel.rad > (6371-bot)*1000);
ivsv = lmodel.vsv(radind);
ivsh = lmodel.vsh(radind);
ivpv = lmodel.vpv(radind);
ivph = lmodel.vph(radind);
irad = lmodel.rad(radind);
irho = lmodel.rho(radind);
vsz = (6371-irad/1000);

depthlim = [0 45];
figure(1)
clf
subplot(1,3,1)
hold on

if islake
    colr = 'b';
else
    colr = 'r';
end

plot(ivpv,vsz,'-b','linewidth',2)
plot(ivsv,vsz,'-r','linewidth',2)
plot(ivsv,vsz,'r.','markersize',25)
plot(ivpv,vsz,'b.','markersize',25)
plot(irho,vsz,'-g','linewidth',2)
plot(irho,vsz,'g.','markersize',25)

set(gca,'ydir','reverse')
xlabel('V(m/s)')
ylabel('Depth (km)')
ylim(depthlim)
xlim([0 8000])
grid on
box on
frechfile = [cardid,'_frech.mat'];

if isrestart
    system(['rm ',frechfile]);
end

if ~exist(frechfile,'file')

frun_mineos([cardid,'.card'],1);
[frech_T,frech_S] = fmk_kernels_linear(cardid,runpath,periods,0);

if ~isempty(frech_S(1).rad)
save(frechfile,'frech_S','lmodel')
end
else
    load(frechfile)
end


%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results for ambient noise
%%%%%%%%%%%%%%%%%%%%%%%%
xaxes = [0e-7 4e-7];
% if islake
% %     CC = cptcmap('GMT_temperature','ncol',length(ANperiods));
% % CC = cbrewer('div','RdYlBu',length(periods));
% CC = jet(length(periods));
% else
% %     CC = cptcmap('GMT_split','ncol',length(periods));
%     CC = cbrewer('div','RdYlGn',length(periods));
% % CC = cptcmap('GMT_temperature','ncol',length(ANperiods));
% end

CC = flipud(cbrewer('seq','Blues',length(ANperiods)+2));
fz = (6371-frech_S(1).rad/1000);
figure(1)
subplot(1,3,2)
hold on
grid on
box on
% plot([0 0],[0 150],'--','color','k')
for ip = 1:length(ANperiods)
   clr = CC(ip,:);
plot(frech_S(ip).vsv,fz,'-','linewidth',2,'color',clr)
lgnd{ip} = [num2str(frech_S(ip).per),' s'];
end
legend(lgnd,'location','southeast')
set(gca,'ydir','reverse')
ylim(depthlim)
xlim(xaxes)

%%%%%%%%%%%%%%%%%%%%%%%%
% Plot all periods
%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if islake
%     CC = cbrewer('div','RdYlBu',length(periods));
% % CC = cptcmap('GMT_temperature','ncol',length(periods));
% else
% %     CC = cptcmap('GMT_split','ncol',length(periods));
%     CC = cbrewer('div','RdYlGn',length(periods));
% % CC = cptcmap('GMT_temperature','ncol',length(periods));
% end
CC = cbrewer('seq','Reds',length(periods)-length(ANperiods)+1);

fz = (6371-frech_S(1).rad/1000);
figure(1)
subplot(1,3,3)
hold on
grid on
box on
count = 0;
clear lgnd
for ip = length(ANperiods)+1:length(periods)
    count = count+1;
   clr = CC(count+1,:);
plot(frech_S(ip).vsv,fz,'-','linewidth',2,'color',clr)
lgnd{count} = [num2str(frech_S(ip).per),' s'];
% lgnd{ip} = [num2str(frech_S(ip).per),' s'];

end
set(gca,'ydir','reverse')
ylim([0 150])
% xlim(xaxes)
xlim([0e-8 3e-8])
legend(lgnd,'location','southeast')



isprint = 1;
if isprint
           set(gcf,'PaperPositionMode','manual');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperPosition',[-0.1 .025 12 8]);
print(figure(1),'-dpsc',['Sensitivity_',cardid,'.ps'])
end
