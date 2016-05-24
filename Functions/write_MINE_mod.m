% write_MINE_mod
% Takes in an input file and writes out finely sampled model card for
% MIENOS
% We only care about the top 250 km because the rest of the model will be
% identical to the anisotropic PREM model we already have
%
% NJA, 2014]
%
% BEWARE -- modifications were made to simplify things for a linear test


function [] = write_MINE_mod(bot,sv_sl,sh_sl,sz,pv_sl,ph_sl,CARD)
% clear
% Controls
isfigure = 0;
yaxis=[0 400];
% turn off warnings momentarily because they are annoying!
warning('off','all');

setup_parameters;
rho = param.rho;
% 
% sed_rho = rho(1);
% crust_rho = rho(2);
% mant_rho = rho(3);
% 

%% Initialize Arrays
% qkap_nsl = ones(length(sz),1);
% qmu_nsl = ones(length(sz),1);
% rho_nsl = ones(length(sz),1);
% eta_nsl = ones(length(sz),1);

%% Get information from original model card
% modelcard = 'EAR_mod.mat';
% load(modelcard)

card = read_model_card('goc.card');

%% Find the depths greater than which we will keep everything the same
dist = 100; % distance overwhich the transition from perturbed velocities to PREM velocities
dp = find(card.z > bot+dist); % dp = deep -- used to be dp = find(card.rad > bot+1);
z_dp = card.z(dp);
rad_dp = card.rad(dp);
sv_dp = card.vsv(dp);
sh_dp = card.vsh(dp);
pv_dp = card.vpv(dp);
ph_dp = card.vph(dp);
% eta_dp = card.eta(dp);
% qkap_dp = card.qkap(dp);
% rho_dp = card.rho(dp);
% qmu_dp = card.qmu(dp);

if isfigure
% figure(19)
% clf
% subplot(2,2,1)
% 
% plot(sv_dp,z_dp,'-k');
% subplot(2,2,2)
% plot(sh_dp,z_dp,'-b');
% subplot(2,2,3)
% plot(pv_dp,z_dp,'-r');
% subplot(2,2,4)
% plot(ph_dp,z_dp,'-m');
end

% Fix NaN in ph
dum = find(isnan(ph_dp)==1);
ph_dp(dum) = 0;

% %% Find the depths for which we are changing things
% sl = find(card.z <= bot+dist);
% z_sl = card.z(sl);

%% Smooth the transition from perturbed velocities to PREM
% Want to create a shallow linear gradient between where the splines stop
% and PREM starts.

ndp = find(card.z > bot & card.z <=bot+dist);

z_tr = card.z(ndp);
rad_tr = card.rad(ndp);
% qmu_tr = card.qmu(ndp);
% rho_tr = card.rho(ndp);
% qkap_tr = card.qkap(ndp);
% eta_tr = card.eta(ndp);
sv_tr = card.vsv(ndp);
sh_tr = card.vsh(ndp);
pv_tr = card.vpv(ndp);
ph_tr = card.vph(ndp);

clear P
%sv
P = polyfit([sz(1) sz(1)+dist],[sv_sl(1) sv_tr(1)],1);
y = polyval(P,z_tr);
sv_ntr = y;

clear P
%sh
P = polyfit([sz(1) sz(1)+dist],[sh_sl(1) sh_tr(1)],1);
y = polyval(P,z_tr);
sh_ntr = y;

clear P
%pv
P = polyfit([sz(1) sz(1)+dist],[pv_sl(1) pv_tr(1)],1);
y = polyval(P,z_tr);
pv_ntr = y;

clear P
%ph
P = polyfit([sz(1) sz(1)+dist],[ph_sl(1) ph_tr(1)],1);
y = polyval(P,z_tr);
ph_ntr = y;

% % rho
% P = polyfit([sz(end) sz(end)+dist],[rho_sl(end) rho_tr(1)],1);
% y = polyval(P,z_tr);
% rho_ntr = y;




if isfigure
    
    figure(22)
    clf
    subplot(1,2,1)
    hold on
    plot(sv_sl,sz','-r','linewidth',2);
    plot(sv_dp,z_dp,'-b','linewidth',2);
    plot(sv_ntr,z_tr,'ok','linewidth',2);
    plot([sv_sl(1) sv_tr(1)],[sz(1) sz(1)+dist],'om','linewidth',2)
    ylim(yaxis)
    set(gca,'ydir','reverse')
    
    subplot(1,2,2)
    hold on
    plot(sh_sl,sz','-r','linewidth',2);
    plot(sh_dp,z_dp,'-b','linewidth',2);
    plot(sh_ntr,z_tr,'-k','linewidth',2);
    ylim(yaxis)
    set(gca,'ydir','reverse')
    
    figure(23)
    clf
    subplot(1,2,1)
    hold on
    plot(pv_sl,sz','-g','linewidth',2);
    plot(pv_dp,z_dp,'-b','linewidth',2);
    plot(pv_ntr,z_tr,'ok','linewidth',2);
    ylim(yaxis)
    set(gca,'ydir','reverse')
    
    subplot(1,2,2)
    hold on
    plot(ph_sl,sz','-g','linewidth',2);
    plot(ph_dp,z_dp,'-b','linewidth',2);
    plot(ph_ntr,z_tr,'ok','linewidth',2);
    ylim(yaxis)
    set(gca,'ydir','reverse')
    
end

nsv = [sv_dp',sv_ntr',sv_sl'];
nsh = [sh_dp',sh_ntr',sh_sl'];
npv = [pv_dp',pv_ntr',pv_sl'];
nph = [ph_dp',ph_ntr',ph_sl'];

nsz = [z_dp',z_tr',sz'];
if isfigure
figure(17)
clf
hold on
plot(nsv,nsz,'-k','linewidth',2)
plot(nsh,nsz,'-r','linewidth',2)
end
% 
% size(nsv)
% size(card.rho)
% return
% Now find all discontinuities - dont need to worry about this happening
% here but you can check if youd like
% [C,ia,ic] = unique(z_dp);
% disc_dp = setdiff((1:length(z_dp)),ia);
% %
% 
% % Shallow depths of old model
% % [C,ia,ic] = unique(z_sl);
% % disc_sl = setdiff((1:length(z_sl)),ia);
% 
% % Shallow depths of new model
% [C,ia,ic] = unique(sz);
% disc_sl = setdiff((1:length(sz)),ia);
% 
% z_disc_dp = z_dp(disc_dp)';
% rad_disc_dp = (6371-z_disc_dp)*1000;
% z_disc_sl = sz(disc_sl)';
% rad_disc_sl = (6371-z_disc_sl)*1000;

% disp(rad_disc_sl);
% if length(disc_sl) < 1
% %     disp(disc_dp)
%     error('No shallow disctontinuities found!')
% end
% 


%
% rad_disc = [rad_disc_dp,rad_disc_sl];
%


% Convert depths to radius
% rad_sl = (6371-sz)*1000;

%% Now begin to write things out to the new model card

fid=fopen(CARD,'w');

% First the header information
% L1 - Model Card Name
fprintf(fid,'%s\n',CARD);

% L2 - ifanis, tref, ifdeck
ifanis=1;
trec = -1;
ifdeck = 1;

fprintf(fid, '%i\t%f\t%i\n',[ifanis trec ifdeck]);

% L3 - N, nic, noc
% N = length(sz)+length(z_disc_all)+length(z_dp);

% N = length(rad_dp)+length(rad_sl)+length(rad_tr);
N = length(card.rho);
nic = 63; % Get from model card -- inner core
noc = 177; % Get from model card -- outer core

fprintf(fid,'%3i\t%2i\t%2i\n',[N nic noc]);

% L4 until end - r, rho, vpv, vsv, qkapp, qshear, vph, vsh, eta
% Remember that r means radius (not depth!)
% Everything else is in m/s (not km/s!)

count = 0;

% First the deeper (unchanged) parts of the model
for id = 1:length(card.rho)
    del = ' ';
    fprintf(fid,'%8.0f%9.2f%9.2f%9.2f%9.1f%9.1f%9.2f%9.2f%9.5f\n'...
        ,[card.rad(id) card.rho(id) npv(id) nsv(id) card.qkap(id) card.qmu(id) nph(id) nsh(id) card.eta(id)]);
    count = count+1;
end

% Now the transition from PREM to our model

% for itr = 1:length(rad_tr)
%     del = ' ';
%     fprintf(fid,'%8.0f%9.2f%9.2f%9.2f%9.1f%9.1f%9.2f%9.2f%9.5f\n'...
%         ,[rad_tr(itr) rho_ntr(itr) pv_ntr(itr) sv_ntr(itr) qkap_tr(itr) qmu_tr(itr) ph_ntr(itr) sh_ntr(itr) eta_tr(itr)]);
%     count = count+1;
% end
% 
% % Then the shallow portion -- this is going to be difficult to figure out
% % how to handle the discontinuities
% for isl = 1:length(rad_sl)
%     is = length(rad_sl)+1-isl;
%     
%     fprintf(fid,'%8.0f%9.2f%9.2f%9.2f%9.1f%9.1f%9.2f%9.2f%9.5f\n'...
%         ,[rad_sl(is) rho_sl(is) pv_sl(is) sv_sl(is) qkap_nsl(is) qmu_nsl(is) ph_sl(is) sh_sl(is) eta_nsl(is)]);
%     count = count+1;
% end

fclose(fid);

if count ~= N
    disp(['N : ',num2str(N),' COUNT : ',num2str(count)]);
    error('Mismatch in layer count!');
end

%turn the warnings back on b/c they can be useful in some cases
warning('on','all')
