% Inversion code for radial anisotropy modified from a bunch of messy messy
% code I wrote earlier in the year.
%
% NJA, November 22, 2014
%

clear

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

icardid = 'goc';

% Parameters for inversion
% maxiter = param.maxiter;
maxiter = 5
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

% Load in Observed Data
% Phase velocity
% TFILE = ['./data/helmholtz_stack_',component,'T.mat'];
% SFILE = ['./data/helmholtz_stack_',component,'Z.mat'];
% load(TFILE);
% Lavg = avgphv;
% clear avgphv;
% load(SFILE);
% Ravg = avgphv;


% Load phase velocity information
ddata = param.data;
load(ddata);

%% Start Inversion
ilat = 5;
ilon = 5;

count = 0;

% Get initial model from goc card

lmodel = read_model_card([icardid,'.card']);

radind = find(lmodel.rad > (6371-bot)*1000);
ivsv = lmodel.vsv(radind);
ivsh = lmodel.vsh(radind);
ivpv = lmodel.vpv(radind);
ivph = lmodel.vph(radind);
irad = lmodel.rad(radind);



while ilat <= Nx
    ilat = ilat+1;
    while ilon <= Ny
        
        % initialize important vectors
        iter = 1;
        count = count+1;
        ilon = ilon+1;
        
        cvsv = ivsv;
        cvsh = ivsh;
        cvpv = ivpv;
        cvph = ivph;
        
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
            continue
        end

        
        
        %% Create kernels
        % first run mineos once
%         frun_mineos([icardid,'.card']);
        
%         return
        CC = jet(maxiter);
        figure(77)
        clf
        figure(78)
        clf
        subplot(1,2,1)
        plot(ivsv,irad,'-k','linewidth',2)
        subplot(1,2,2)
        plot(ivsh,irad,'-k','linewidth',2)
        

        while iter <= maxiter
            %% Create Initial Model

            if iter == 1
                cardid = icardid;
                frun_mineos([cardid,'.card'],1);
                [frech_T,frech_S] = fmk_kernels_linear(icardid,runpath,periods);
                
                figure(73)
                clf
                hold on
                plot(periods,obs.Lphv,'or','linewidth',2)
                plot(periods,obs.Rphv,'ob','linewidth',2)
%                 pause
            else
                cardid = [icardid,num2str(iter)];
                write_MINE_mod(bot,cvsv,cvsh,6371-(irad/1000),cvpv,cvph,[cardid,'.card']);

            end
            
            %% ----- Run Mineos ------
            tic
            frun_mineos([cardid,'.card'],0);
            toc
            
            % Compare q correction and non-q corrected.
            
            
            
            % Grab the estimated phase velocities
            clear TPHV TCPER SPHV SCPER

            TYPE = 'T';
            [TCPER,TPHV] = calc_fundphv(periods,TYPE);
%             [TCPER,TPHV,TUPER,TGRV,TMODE] = calc_fundCU(Cperiods,Uperiods,TYPE);
            TYPE = 'S';
            [SCPER,SPHV] = calc_fundphv(periods,TYPE);
%             [SCPER,SPHV,SUPER,SGRV,SMODE] = calc_fundCU(Cperiods,Uperiods,TYPE);
            
            forward.tphv = TPHV;
            forward.tper = TCPER;
            forward.sphv = SPHV;
            forward.sper = SCPER;
            
            %% Run the inversion
%             save('mest_test.mat','obs.Lphv','obs.Lstd','obs.Rphv','obs.Rstd','forward','frech_T','frech_S');
            mest = finverse_aniso_linear(obs.Lphv,obs.Lstd,obs.Rphv,obs.Rstd,forward,frech_T,frech_S,periods);
            
            %% update the model parameters
            % First SV
            pertSV = mest.SV;
             vsv = cvsv+pertSV;
            
            % Next SH
            pertSH = mest.SH;
            vsh = cvsh+pertSH;
            

            % Perturb the p wave 
            vpv = cvpv.*(1+pertSV./vsv);
            vph = cvph.*(1+pertSH./vsh);
%             coef_pvn = coef_pv.*(1+pertSV'./coef_sv);
%             coef_phn = coef_ph.*(1+pertSH'./coef_sh);
            
            % update vectors for next run
            cvsv = vsv;
            cvsh = vsh;
            cvpv = vpv;
            cvph = vph;
            
            figure(77)
            subplot(1,2,1)
            hold on
            plot(forward.tper,forward.tphv,'o','markerfacecolor',CC(iter,:),'markersize',7,'markeredgecolor','k')
            subplot(1,2,2)
            hold on
            plot(forward.sper,forward.sphv,'o','markerfacecolor',CC(iter,:),'markersize',7,'markeredgecolor','k')

            figure(78)
            subplot(1,2,1)
            hold on
            plot(vsv,irad,'-k','color',CC(iter,:))
            xlim([2700 4700]);
            ylim([6070000 6371000])
            title('Vsv');
            subplot(1,2,2)
            hold on
            plot(vsh,irad,'-k','color',CC(iter,:));
            xlim([3000 4800])
            ylim([6070000 6371000])
            title('Vsh')
                        iter = iter+1;
%                         pause
        end % max iter
        figure(77)
        subplot(1,2,1)
        hold on
        plot(periods,obs.Lphv,'or','linewidth',2,'markerfacecolor','k')
        xlim([5 100])
        title('Love Wave')
        subplot(1,2,2)
        hold on
        plot(periods,obs.Rphv,'ob','linewidth',2,'markerfacecolor','k')
        title('Rayleigh Wave')
        xlim([5 100])
        %         result.obs.Lphv = obs.Lphv;
        %         result.obs.Lstd = obs.Lstd;
        %         result.obs.Rphv = obs.Rphv;
        %         result.obs.Rstd = obs.Rstd;
        % plot the results of hte inversion
%         return

        print(figure(77),'-dpsc',['Phv_',num2str(xi(ilat,ilon)),'x',num2str(yi(ilat,ilon)),'.ps']);
        print(figure(78),'-dpsc',['Vs_',num2str(xi(ilat,ilon)),'x',num2str(yi(ilat,ilon)),'.ps']);
        pause
    end % ilon
end % ilat
