% Inversion code to invert phase velocities for perturbations to spline
% coefficients
%
% NJA, 11/22/2014
%
% Based off of code from 5/2014
%

function mest = finverse_aniso_linear(Lphv,Lstd,Rphv,Rstd,forward,frech_T,frech_S,periods)
% Lphv = obs.Lphv;
% Rphv = obs.Rphv;
% Lstd = obs.Lstd;
% Rstd = obs.Rstd;
% % 

isfigure = 0;

setup_parameters;
% periods = param.Cperiods;
bot = param.bot;
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
    
%     dCL(ip) = (tphv(tind)-Lphv(ip));
%     dCR(ip) = (sphv(sind)-Rphv(ip));
      dCL(ip) = (Lphv(ip)-tphv(tind));
      dCR(ip) = (Rphv(ip)-sphv(tind));
end

dobs = [dCL,dCR]';

nanind = find(isnan(dobs)==1);

dobs(nanind) = 0;

% change to m/s
dobs = dobs*1000;

% disp(dobs)

%% Build G matrix - pull in coefficient specific sensitivities from the
% perturb_mod script
% This will be organized with the Love wave sensitivities first and then
% the VSV sensitivities
% 
% dTSVdCOEF = kern.TSV;
% dTSHdCOEF = kern.TSH;
% dSSVdCOEF = kern.SSV;
% dSSHdCOEF = kern.SSH;
% 
% % check to make sure matrices look correct
% if size(dTSVdCOEF) ~= size(dSSVdCOEF)
%     disp('Error! Kernel Matrices dont match!')
%     error;
% elseif size(dTSHdCOEF) ~= size(dSSHdCOEF)
%         disp('Error! Kernel Matrices dont match!')
%     error;
% end
% 
% [Ncoef Nper] = size(dTSHdCOEF);

% Find only those sensitiviites ar at depth we care about.  

tindrad = find(frech_T(1).rad > (6371-bot)*1000);
sindrad = find(frech_S(1).rad > (6371-bot)*1000);

% length(tindrad)
% length(sindrad)
% % 

GG = zeros(ip*2,length(tindrad)*2);

% size(GG)

% Sensitivities from Love
for ip = 1:length(frech_T)
    GG(ip,1:length(tindrad)) = frech_T(ip).vsv(tindrad)';
    GG(ip,length(tindrad)+1:length(tindrad)*2) = frech_T(ip).vsh(tindrad);
end

% Sensitivities from Rayleigh
for ip = 1:length(frech_S)
    GG(ip+length(frech_T),1:length(sindrad)) = frech_S(ip).vsv(sindrad);
%     GG(ip+length(frech_T),length(sindrad)+1:length(sindrad)*2) = zeros(length(sindrad),1);
%     GSH(ip+length(frech_T),:) = frech_S(ip).vsh(sindrad);
    
end


% Make the smoothness matrix
[M N] = size(GG);

% GG = GG/(1E6);

Wm = zeros(N,N);

D = zeros(N,N);
for in = 2:N-1
    D(in,in-1) = 1;
    D(in,in) = -2;
    D(in,in+1) =1;
end
D(1,1) = -1;
D(1,2) = 1;
D(N,N-1) = -1;
D(N,N) = 1;
% N
% M
% for in = 1:(N-2)
% %         if in == N
% %             Wm(N,N-1) = -1;
% %             Wm(N,N) = 1;
% %         else
% if in == 1 || in == 2
% %     Wm(in,in) = 1;
% %     Wm(in,in+1) = -1;
% else
%             Wm(in,in) = 1;
%             Wm(in,in+1) = -2;
%             Wm(in,in+2) = 1;
% end
% end

for in = 1:N-1
    Wm(in,in) = 1;
    Wm(in,in+1) = -1;
end

% Wm
% Wm(N-1,N-2) = 1;
% Wm(N-1,N-1) = -1;
% Wm(N,N-1) = 1;
% Wm(N,N) = -1;
% Make the error weight matrix

We = zeros(ip*2,ip*2);

for ip = 1:length(periods)
    
    if Lstd(ip) == 0
        LTstd = 0.15;
    else
        LTstd = Lstd(ip);
    end
    if Rstd(ip) == 0
        RTstd = 0.15;
    else
        RTstd = Rstd(ip);
    end
    We(ip,ip) = LTstd;
    We(ip+length(periods),ip+length(periods)) = RTstd;
end

We = We;
% We

% return

if isfigure
% double check the G matrix via observation

figure(1)
clf

CC = jet(length(frech_T));
for ip=1:length(frech_T)
subplot(1,2,1)
hold on
    plot(GG(ip+length(frech_T),1:length(sindrad)),frech_S(ip).rad(sindrad),'-k','color',CC(ip,:));
    plot(GG(ip,1:length(tindrad)),frech_T(ip).rad(tindrad),'--','color',CC(ip,:))
   ylim([(6371-bot)*1000 6371*1000])
   title('Vsv')
   xlim([0 3E-8])
    subplot(1,2,2)
hold on
plot(GG(ip,length(tindrad)+1:length(tindrad)*2),frech_T(1).rad(tindrad),'--','color',CC(ip,:));
% plot(GSH(ip+length(frech_T),:),frech_S(ip).rad(sindrad),'-','color',CC(ip,:));
ylim([(6371-bot)*1000 6371*1000])
title('Vsh')
end


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

%% Genearte the weight matrix -- based off of ratios of standard deviation
% Most accurate measurements will have 

% Inverse the matrix to get mest (perturbations to the coeffi
% mest.SV = (GSV'*GSV)\(GSV'*dobs);

% mest.SV = lsqr(GSV,dobs);
% mest_all = lsqr(GG,dobs);

% weighted damped least squares:

mhat = zeros(length(tindrad)*2,1);

epsilon = 1E-13;
% % 
% % size(GG)
% % size(We)
% 
% % size(GG'*We*GG)
% % size(eta^2*Wm)
% mest_all1 = [GG'*We*GG+eta^2*Wm];
% mest_all2 = [GG'*We*dobs+eta^2*Wm*mhat];
% mest_all = mest_all1\mest_all2;
% mest_all = mest_all/(1e6);
% mest.SV = mest_all(1:length(tindrad));
% mest.SH = mest_all(length(tindrad)+1:length(tindrad)*2);

% let's try doing it a different way with the F'Fmest=F'f

% F = [sqrt(We)*GG;eta*D];
% f = [sqrt(We)*dobs;eta*D*mhat];
% 
% % mest = bicg(@weightedleastsquaresfcn,F'*f,tol,maxit);
% mest_all = (F'*F)\(F'*f);
% mest_all = mest_all/1E9;
% disp(mest_all)

% lets try just the weighted least squares
clear mest_all
% mest_all = (GG'*We*GG)\(GG'*We*dobs);

% lets try with minimum length
clear mest_all
GMG = GG'/(GG*GG'+epsilon*eye(M,M));
mest_all = GMG*dobs;
% mest.SV = mest.SV/(1000^3);
% mest.SV = mest.SV/(10^8);
mest_all = mest_all/(1E6);
mest.SV = mest_all(1:length(tindrad));
mest.SH = mest_all(length(tindrad)+1:length(tindrad)*2);
disp('SV')
disp(mest.SV);
disp('SH')
disp(mest.SH);

% mest.SH = (GSH'*GSH)\(GSH'*dobs);
% mest.SH = lsqr(GSH,dobs);
% mest.SH = lsqr(GSH,dCL');
% 
% mest.SH = mest.SH/(1e6);

% mest.SH = mest.SH/(10^8);
% mest.SH = mest.SH/(1000^3);

% disp('SH')
% disp(mest.SH);
% 
% save(savefile,'mest');

%     end % Ny
% end % Nx
