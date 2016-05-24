% Function to read in the initial model card and plot it to double check
% that everything looks ok
%
% varargin used to take in name of model card for iteration step.
% if no name is given then the model is named init_model
% NJA, 2014

function card=read_model_card(CARD)

% datapath = 'data/';

setup_parameters;

bot = param.bot;

isfigure = 0;
% Parameters spcefic to the format of these files
hlines1 = 3;

% Read model card information
% cardname = [datapath,CARD];

% CARD = [CARD,'.card'];
fid = fopen(CARD,'r');
A=textscan...
    (fid,'%f %f %f %f %f %f %f %f %f','headerlines',hlines1);

card.fname=CARD;
card.z = 6371-A{1}/1000;
card.rad = A{1};
card.rho = A{2};
card.vpv = A{3};
card.vph = A{7};
card.vsv = A{4};
card.vsh = A{8};
card.eta = A{9};
card.qmu = A{6};
card.qkap = A{5};

fclose(fid);

% Choose which card you want to look at

if isfigure
yaxis=[0 350];

figure(1)
clf
plot(card.vsv,card.z,'-','color','r','linewidth',2);
hold on
plot(card.vsh,card.z,'-k','linewidth',2);

ylim(yaxis)
xlim([3000 5000])
set(gca,'ydir','reverse','fontsize',16)

hleg = legend('SV','SH');



figure(2)
clf
subplot(2,3,1)
plot(card.vpv,card.z,':','color','r','linewidth',2);
hold on
plot(card.vph,card.z,'-r','linewidth',2);
ylim(yaxis)
set(gca,'ydir','reverse','fontsize',16)
title('Vpv and Vph')
subplot(2,3,2)
plot(card.eta,card.z,':','color','r','linewidth',2);
ylim(yaxis)
set(gca,'ydir','reverse','fontsize',16)
title('Eta')
subplot(2,3,3)
plot(card.qmu,card.z,':','color','r','linewidth',2);
ylim(yaxis)
set(gca,'ydir','reverse','fontsize',16)
title('Qmu')
subplot(2,3,4)
plot(card.qkap,card.z,':','color','r','linewidth',2);
ylim(yaxis)
set(gca,'ydir','reverse','fontsize',16)
title('Q Kappa')
subplot(2,3,5)
plot(card.rho,card.z,':','color','r','linewidth',2);
ylim(yaxis)
set(gca,'ydir','reverse','fontsize',16)
title('Density')


% Let's plot the anisotropy (SH^2/SV^2)

vsv2 = card.vsv.^2;
vsh2 = card.vsh.^2;

psi = (vsh2./vsv2);

figure(3)
clf
ydum = 0:10:600;
xdum = ones(size(ydum));
plot(psi,card.rad,'-b','linewidth',2)
hold on
plot(xdum,ydum,'.g','markersize',10);
ylim(yaxis)
xlim([0.95 1.1])
set(gca,'ydir','reverse','fontsize',16)
title('Anisotropy');
end
% figure(2)
% clf
% subplot(1,2,1)
% plot(card(3).vsv,card(3).rad,':','color','r','linewidth',2);
% hold on
% plot(card(3).vsh,card(3).rad,'-r','linewidth',2);
% ylim([0 100])
% set(gca,'ydir','reverse','fontsize',16)
% subplot(1,2,2)
% plot(mod(2).t,mod(2).grv,':','color','r','linewidth',2);
% hold on
% plot(mod(2).t,mod(2).phv,'-','color','r','linewidth',2);
% xlim([0 150])
% set(gca,'fontsize',16);