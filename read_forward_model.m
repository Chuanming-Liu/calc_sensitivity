% Program to plot MINEOS model cards
% NJA, 2014
%
% 10/27/2014 Modified to handle either longer period S models (S0_50) or
% shorter period S models (S50_150).

function mode = read_forward_model(CARD,TYPE)
% tic
% CARD='EARS.card';
% TYPE = 'T';
% FREQ = '';

isfigure = 0;

% Get useful information for parameter file
setup_parameters;

DATAPATH = param.DATAPATH;
SID = param.SID;
TID = param.TID;
% dum = strsplit(CARD,'.');
% 
% CARD = dum{1};


if strcmp(TYPE,'S') == 1
    
    CASC = [DATAPATH,CARD,'.',SID,'.asc'];

elseif strcmp(TYPE,'T') == 1
    CASC = [DATAPATH,CARD,'.',TID,'.asc'];
else
    disp('Type does not exist! Use T or S')
end

% Parameters spcefic to the format of these files
hlines1 = 5;
hlines2 = 6;

% Read model card information
cardname = CASC;
fid = fopen(cardname,'r');
A=textscan...
    (fid,'%f %f %f %f %f %f %f %f %f %f','headerlines',hlines1);

card.fname=cardname;
card.lev = A{1};
card.rad = 6371-A{2}/1000;
card.rho = A{3};
card.vpv = A{4};
card.vph = A{5};
card.vsv = A{6};
card.vsh = A{7};
card.eta = A{8};
card.qmu = A{9};
card.qkap = A{10};


% Pick up where we left off but with a different format
B=textscan(fid,'%d %s %d %f %f %f %f %f %f','headerlines',hlines2);

mode.fname=cardname;
mode.n = B{1};
mode.l = double(B{3});
mode.wrad = B{4};
mode.w = B{5};
mode.t = B{6};
mode.grv = B{7};
mode.q = B{8};
mode.rayq = B{9};

fclose(fid);

for im=1:length(mode.n)
mode.phv(im) = mode.w(im)*2*pi*(6371)/sqrt(mode.l(im)*(mode.l(im)+1))/1000;
end

savefile = [CARD,'.',TYPE,'.mat'];

save(savefile,'card','mode');
% toc
% Choose which card you want to look at

yaxis=[0 500];

% figure(1)
% clf
% subplot(1,2,1)
% plot(card.vsv,card.rad,':','color','k','linewidth',2);
% hold on
% plot(card.vsh,card.rad,'-k','linewidth',2);
% ylabel('Depth');
% xlabel('Vs');
% ylim(yaxis)
% set(gca,'ydir','reverse','fontsize',16)
% title('Vsv vs. Vsh');
% subplot(1,2,2)
% plot(mod.t,mod.grv,':','color','k','linewidth',2);
% xlabel('Period (s)')
% ylabel('Group Velocity');
% hold on
% % plot(mod.t,mod.phv,'-','color','k','linewidth',2);
% xlim([0 150])
% set(gca,'fontsize',16);
% title('Group Velocity against Period');

% Find all fundamental modes

fundind = find(mode.n == 0);

fund_grv = mode.grv(fundind);
fund_t = mode.t(fundind);
fund_phv = mode.phv(fundind);

fund.grv = fund_grv;
fund.per = fund_t;
fund.phv = fund_phv;

if isfigure
figure(3)
clf
plot(fund_t,fund_grv,'-b','linewidth',2)
hold on
plot(fund_t,fund_phv,'--b','linewidth',2)
xlim([0 300])
set(gca,'fontsize',16);
title(['Fundamental Mode: ',FREQ]);




figure(2)
clf
subplot(2,3,1)
plot(card.vpv,card.rad,':','color','r','linewidth',2);
hold on
plot(card.vph,card.rad,'-r','linewidth',2);
ylim(yaxis)
set(gca,'ydir','reverse','fontsize',16)
title('Vpv and Vph')
subplot(2,3,2)
plot(card.eta,card.rad,':','color','r','linewidth',2);
ylim(yaxis)
set(gca,'ydir','reverse','fontsize',16)
title('Eta')
subplot(2,3,3)
plot(card.qmu,card.rad,':','color','r','linewidth',2);
ylim(yaxis)
set(gca,'ydir','reverse','fontsize',16)
title('Qmu')
subplot(2,3,4)
plot(card.qkap,card.rad,':','color','r','linewidth',2);
ylim(yaxis)
set(gca,'ydir','reverse','fontsize',16)
title('Q Kappa')
subplot(2,3,5)
plot(card.rho,card.rad,':','color','r','linewidth',2);
ylim(yaxis)
set(gca,'ydir','reverse','fontsize',16)
title('Density')
% pause
end