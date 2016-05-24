% Function based off of the write_MINE_mod function to write out Mineos
% card files but modified to handle discontinuities in the shallowest parts
% 
% NJA, 4/20/2016


function [] = create_Mineos_cards(sv_sl,zz,pv_sl,qk_sl,qs_sl,rho_sl,CARD)
% clear
% Controls
isfigure = 0;
yaxis=[0 400];

% turn off warnings momentarily because they are annoying!
warning('off','all');

setup_parameters;
rho = param.rho;

sh_sl = sv_sl;
ph_sl = pv_sl;

bot = zz(1);
rad_sl = (6371-zz)*1000;
%% Get information from original model card

card = read_model_card('goc.card');

%% Find the depths greater than which we will keep everything the same
dist = 50; % distance overwhich the transition from perturbed velocities to PREM velocities
dp = find(card.z > bot+dist); % dp = deep -- used to be dp = find(card.rad > bot+1);
z_dp = card.z(dp);
rad_dp = card.rad(dp);
sv_dp = card.vsv(dp);
pv_dp = card.vpv(dp);
eta_dp = card.eta(dp);
qk_dp = card.qkap(dp);
rho_dp = card.rho(dp);
qs_dp = card.qmu(dp);

if isfigure
figure(19)
clf
subplot(1,2,1)
plot(sv_dp,z_dp,'-r');
subplot(1,2,2)
plot(pv_dp,z_dp,'-b');
end


%% Smooth the transition from perturbed velocities to PREM
% Want to create a shallow linear gradient between where the splines stop
% and PREM starts.

ndp = find(card.z > bot & card.z <=bot+dist);

z_tr = card.z(ndp);
rad_tr = card.rad(ndp);
qs_tr = card.qmu(ndp);
rho_tr = card.rho(ndp);
qk_tr = card.qkap(ndp);
eta_tr = card.eta(ndp);
sv_tr = card.vsv(ndp);
pv_tr = card.vpv(ndp);

clear P
ind = 1;
%sv
P = polyfit([zz(ind) zz(ind)+dist],[sv_sl(ind) sv_tr(ind)],1);
y = polyval(P,z_tr);
sv_ntr = y;
sh_ntr = sv_ntr;

clear P
%pv
P = polyfit([zz(ind) zz(ind)+dist],[pv_sl(ind) pv_tr(ind)],1);
y = polyval(P,z_tr);
pv_ntr = y;
ph_ntr = pv_ntr;

if isfigure
    
    figure(22)
    clf
    subplot(1,2,2)
    hold on
    plot(sv_sl,zz','-r','linewidth',2);
    plot(sv_dp,z_dp,'-b','linewidth',2);
    plot(sv_ntr,z_tr,'ok','linewidth',2);
    plot([sv_sl(1) sv_tr(1)],[zz(1) zz(1)+dist],'om','linewidth',2)
    ylim(yaxis)
    set(gca,'ydir','reverse')
    subplot(1,2,1)
    hold on
    plot(pv_sl,zz','-g','linewidth',2);
    plot(pv_dp,z_dp,'-b','linewidth',2);
    plot(pv_ntr,z_tr,'ok','linewidth',2);
    ylim(yaxis)
    set(gca,'ydir','reverse')
    
end

nsv = [sv_dp',sv_ntr',sv_sl'];
nsh = nsv;
npv = [pv_dp',pv_ntr',pv_sl'];
nph = npv;
rho = [rho_dp',rho_tr',rho_sl'];
qk = [qk_dp',qk_tr',qk_sl'];
qs = [qs_dp',qs_tr',qs_sl'];
eta = [eta_dp',eta_tr',ones(size(qs_sl'))];
rad = [rad_dp',rad_tr',rad_sl'];

nsz = [z_dp',z_tr',zz'];
if isfigure
figure(17)
clf
hold on
subplot(1,2,1)
plot(npv,nsz,'-b','linewidth',2)
subplot(1,2,2)
plot(nsv,nsz,'-r','linewidth',2)
end


% Convert depths to radius
% rad_sl = (6371-sz)*1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Now begin to write things out to the new model card %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen(CARD,'w');

% First the header information
% L1 - Model Card Name
fprintf(fid,'%s\n',CARD);

% L2 - ifanis, tref, ifdeck
ifanis=0;
trec = -1;
ifdeck = 1;

fprintf(fid, '%i\t%f\t%i\n',[ifanis trec ifdeck]);

% L3 - N, nic, noc
% N = length(sz)+length(z_disc_all)+length(z_dp);

% N = length(rad_dp)+length(rad_sl)+length(rad_tr);
% N = length(card.rho);
N = numel(npv);
nic = 63; % Get from model card -- inner core
noc = 177; % Get from model card -- outer core

fprintf(fid,'%3i\t%2i\t%2i\n',[N nic noc]);

% L4 until end - r, rho, vpv, vsv, qkapp, qshear, vph, vsh, eta
% Remember that r means radius (not depth!)
% Everything else is in m/s (not km/s!)

count = 0;

% First the deeper (unchanged) parts of the model
for id = 1:N
    del = ' ';
    fprintf(fid,'%8.0f%9.2f%9.2f%9.2f%9.1f%9.1f%9.2f%9.2f%9.5f\n'...
        ,[rad(id) rho(id) npv(id) nsv(id) qk(id) qs(id) nph(id) nsh(id) eta(id)]);
    count = count+1;
end

fclose(fid);

if count ~= N
    disp(['N : ',num2str(N),' COUNT : ',num2str(count)]);
    error('Mismatch in layer count!');
end

%turn the warnings back on b/c they can be useful in some cases
warning('on','all')
