% Calculate ratio between individual spline coefficients 
%
% We feed in a model card and then grab the general structure from it
% The velocity ratios will be used to scale the spline coefficients based
% on bulk velocity estiamtes.
% 
% This script is used inside the GS_driver script
%
% NJA, 9/30/2014 - script is based off spline_mod_S
%
% varargin represents the modified spline coefficients to be used when
% iterating. The first variable must for vpv and the second for vph
%
% modified to calculate scales for both the P nad S models

function [Sscale,Pscale,coef_z] = calc_scale(spline,modelmat)

setup_parameters;

bot = param.bot;

isfigure = 0;
yaxis = [0 bot];

% %% Turn on if you want to run it outside of function form
% Calculate the basis functions for the splines
% [spline] = calc_spline(bot,sedh,mohoh);
% modelcard = 'EAR_mod.mat';

warning('off','all')

% Parameters from original model card --- need to change this
omohoh = 19.5;
osedh = 3.5;

% Initialize vectors for the splines
crustv = ones(size(spline(2).b,1),1);
mantv = ones(size(spline(3).b,1),1);

%% Pull velociites gradients the model card
load(modelmat)

rad = card.rad;
ovsv = card.vsv;
ovsh = card.vsh;
ovpv = card.vpv;
ovph = card.vph;

% Create vectors representative of the sediment, crust, and mantle
% velocities

% Sediments
% dum = find(rad < osedh);
% vsv = ovsv;
% vsh = ovsh;
sed_z = mean(0:osedh);
% sed_sv = mean(vsv(dum));
% sed_sh = mean(vsh(dum));

% if sed_sv ~= sed_sh
%     disp('Anisotropy in the sediments!');
% end

Sscale.sed_sv = 1;
Sscale.sed_sh = 1;
Pscale.sed_pv = 1;
Pscale.sed_ph = 1;

%%% Crust %%%%
clear dum vsv vsh trad
dum = find(rad >= osedh & rad <= omohoh);
trad = rad(dum);

% This is going to be messy but I don't care
% Find the beginning and endpoints of the crust
trad = trad(2:end-1);

vsv = ovsv(dum);
vsv = vsv(2:end-1);

vsh = ovsh(dum);
vsh = vsh(2:end-1);

vpv = ovpv(dum);
vpv = vpv(2:end-1);

vph = ovph(dum);
vph = vph(2:end-1);

% Fit a polynomial to the depth and velocities
sP = polyfit(trad,vsv,1);

pP = polyfit(trad,vpv,1);

% Generate the depth intervals we want to look at 
crust_z = trad(end):trad(1);
ind = round(linspace(1,length(crust_z),length(crustv)));
crust_z=crust_z(ind);

% get sv
sy = polyval(sP,crust_z);
crust_z = crust_z';
crust_sv = sy;
crust_sv = crust_sv';

py = polyval(pP,crust_z);
crust_pv = py;
crust_pv = crust_pv';

% Now sh
clear sP sy pP py
sP = polyfit(trad,vsh,1);
sy = polyval(sP,crust_z');
crust_sh = sy;
crust_sh = crust_sh';

pP = polyfit(trad,vph,1);
py = polyval(pP,crust_z');
crust_ph = py;
crust_ph = crust_ph';

% Save scalings to vector
Sscale.crust_sv = [(crust_sv(1)/crust_sv(2)) 1 crust_sv(3)/crust_sv(2)];
Sscale.crust_sh = [(crust_sh(1)/crust_sh(2)) 1 crust_sh(3)/crust_sh(2)];

Pscale.crust_pv = [(crust_pv(1)/crust_pv(2)) 1 crust_pv(3)/crust_pv(2)];
Pscale.crust_ph = [(crust_ph(1)/crust_ph(2)) 1 crust_ph(3)/crust_ph(2)];

% Mantle
clear dum vsv vsh trad vpv vph
dum = find(rad >= omohoh & rad <= (bot+10));
trad = rad(dum);
trad = trad(1:end-1);

vsv = ovsv(dum);
vsv = vsv(1:end-1);
vsh = ovsh(dum);
vsh = vsh(1:end-1);

vpv = ovpv(dum);
vpv = vpv(1:end-1);
vph = ovph(dum);
vph = vph(1:end-1);

% Fit a polynomial to the depth and velocities
sP = polyfit(trad,vsv,4);

pP = polyfit(trad,vpv,4);

% Generate the depth intervals we want to look at
mant_z = trad(end):trad(1);
ind = round(linspace(1,length(mant_z),length(mantv)));
mant_z = mant_z(ind);

% sv
sy = polyval(sP,mant_z);
mant_z = mant_z';
mant_sv = sy;
mant_sv = mant_sv';

py = polyval(pP,mant_z);
mant_pv = py;
mant_pv = mant_pv';

% Now sh
clear sP sy pP py
sP = polyfit(trad,vsh,4);
sy = polyval(sP,mant_z');
mant_sh = sy;
mant_sh = mant_sh';

pP = polyfit(trad,vph,4);
py = polyval(pP,mant_z');
mant_ph = py;
mant_ph = mant_ph';

% Save scalings to a vector

Sscale.mant_sv = [mant_sv(1)/mant_sv(3) mant_sv(2)/mant_sv(3) 1 mant_sv(4)/mant_sv(3) mant_sv(5)/mant_sv(3)];
Sscale.mant_sh = [mant_sh(1)/mant_sh(3) mant_sh(2)/mant_sh(3) 1 mant_sh(4)/mant_sh(3) mant_sh(5)/mant_sh(3)];

Pscale.mant_pv = [mant_pv(1)/mant_pv(3) mant_pv(2)/mant_pv(3) 1 mant_pv(4)/mant_pv(3) mant_pv(5)/mant_pv(3)];
Pscale.mant_ph = [mant_ph(1)/mant_ph(3) mant_ph(2)/mant_ph(3) 1 mant_ph(4)/mant_ph(3) mant_ph(5)/mant_ph(3)];

% Save Depths of Knotts
coef_z = [sed_z;crust_z;mant_z];

% if isfigure
%     figure(67)
%     clf
%     hold on
%     plot(ovsv,rad,'-k','linewidth',2)
%     plot(ovsv,rad,'.k','markersize',12)
%     plot(ovsh,rad,':','color','k','linewidth',2)
%     
%     % Plot the chosen velocity points
%     plot(sed_sh,sed_z,'og','linewidth',2,'markersize',15);
%     plot(mant_sh,mant_z,'ob','linewidth',2,'markersize',15);
%     plot(crust_sh,crust_z,'or','linewidth',2,'markersize',15);
%     
%     plot(sed_sv,sed_z,'xg','linewidth',2,'markersize',15);
%     plot(mant_sv,mant_z,'xb','linewidth',2,'markersize',15);
%     plot(crust_sv,crust_z,'xr','linewidth',2,'markersize',15);
%     ylim(yaxis)
%     set(gca,'ydir','reverse','fontsize',16)
% end



warning('on','all');
