% Inversion code for radial anisotropy modified from a bunch of messy messy
% code I wrote earlier in the year.
%
% NJA, November 22, 2014
%
% WARNING - MODIFIED FOR LINEAR MODELS -- 12/29/2014

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
Uperiods = param.Uperiods;
runpath = param.RUNPATH;
cardid = param.CARDID;

% Parameters for inversion
maxiter = param.maxiter;
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

% Choose parameters for the starting model
sedh1 = sedh(2);
sedv1 = sedv(1);
crustv1 = crustv(2);
mantlev1 = mantlev(2);
mantlepsi1 = mantlepsi(3);
mantstructS1 = mantstructS(1,:);
mantstructP1 = mantstructP(1,:);

% parameters from phase velocity
Cperiods = param.Cperiods;

% turn on use of psi parameterization
ispsi = 1;
ischeckmodel = 0;

maxchi2 = param.maxchi2;

modelmat = 'EAR_mod.mat';

% Set path to executables
set_mineos_path;
% addpath('/home/segment/Natalie/bspline');

% Load in Observed Data
% Phase velocity
TFILE = ['./data/helmholtz_stack_',component,'T.mat'];
SFILE = ['./data/helmholtz_stack_',component,'Z.mat'];
load(TFILE);
Lavg = avgphv;
clear avgphv;
load(SFILE);
Ravg = avgphv;

% Load moho information
load('./data/moho.mat');

%% Start Inversion
ilat = 7;
ilon = 9;

count = 0;
while ilat <= Nx
    ilat = ilat+1;
    while ilon <= Ny
        iter = 1;
        count = count+1;
        ilon = ilon+1;
        
        disp(['PROCESSING ',num2str(xi(ilat,ilon)),'x',num2str(yi(ilat,ilon))]);
        % Get grid point observations
        for ip=1:length(Cperiods) % stupid indexing bc we're skipping 20 s
            obs.Lphv(ip) = Lavg(ip).GV_cor(ilat,ilon);
            obs.Lstd(ip) = Lavg(ip).GV_cor_std(ilat,ilon);
            obs.Rphv(ip) = Ravg(ip).GV_cor(ilat,ilon);
            obs.Rstd(ip) = Ravg(ip).GV_cor_std(ilat,ilon);
        end
%         
%         for ig = 1:length(Uperiods)
%             obs.Rgrv(ig) = double(interp2(MERgrd(ig).lat,MERgrd(ig).lon,MERgrd(ig).grv,xi(ilat,ilon),yi(ilat,ilon)));
%         end
%         obs.RUstd = 0.15*ones(size(obs.Rgrv));
%         
        % check to see if the velocities are empty
        Rnan = find(isnan(obs.Rphv)==1);
        Lnan = find(isnan(obs.Lphv)==1);
%         Gnan = find(isnan(obs.Rgrv)==1);
        
        if length(Rnan) == length(Cperiods)
            disp('No measurements found!')
            continue
        elseif isempty(obs.Rgrv) || isempty(obs.RUstd)
            disp('No group velocity information found!');
            obs.Rgrv = [];
            obs.RUstd = [];
        end
        
        % Get crusth from moho model
        crusth = interp2(moho.xi,moho.yi,moho.depth,xi(ilat,ilon),yi(ilat,ilon));
        crusth = round(crusth);
        ocrusth(count) = crusth;
        disp(['Moho Depth :',num2str(crusth)]);
        if isnan(crusth)
            disp('No moho found! Breaking now ...')
            return;
        end
        
        % Create the splines -- shouldn't change as we iterate
        % Calculate the basis functions for the splines
        [spline] = calc_basis(bot,sedh1,crusth);
        
        % Calculate the scale factors to be used with the new model
        [Sscale,Pscale,coef_z] = calc_scale(spline,modelmat);
        
        %% Check to see if we need to recreate the kernels
        
        if count > 1 && ocrusth(count-1) == crusth
            disp('Previous models kernels will work!')
        else
            % make frechet kernels
            [frech_T,frech_S] = fmk_kernels(cardid,runpath,[Uperiods,Cperiods]);
            
            % make a simple velocity model for the kernels -- needed for
            % conversion
            kcoef.sv = [(Sscale.sed_sv*sedv1),(Sscale.crust_sv*crustv1),(mantstructS1*mantlev1)];
            kcoef.psi = [ones(size(Sscale.sed_sv)),ones(size(Sscale.crust_sv)),(mantstructPSI*mantlepsi1)];
            kcoef.z = coef_z;
            
            for ish = 1:length(kcoef.sv)
                kcoef.sh(ish) = sqrt(kcoef.psi(ish)*(kcoef.sv(ish)^2));
            end
            
            % convet linear model kernels to spline kernels
%             kern = fconv_kern(Cperiods,cardid,frech_T,frech_S,spline,kcoef);
        end
        % return
        % start to iterate
%         return
        CC = spring(maxiter);
        while iter <= maxiter
            %% Run the forward problem
            % completely create the model if this is the first iteration
            if iter == 1
                
                % Create the spline velocity model
                coef_sv = [(Sscale.sed_sv*sedv1),(Sscale.crust_sv*crustv1),(mantstructS1*mantlev1)];
                coef_psi = [ones(size(Sscale.sed_sv)),ones(size(Sscale.crust_sv)),(mantstructPSI*mantlepsi1)];
                for ish = 1:length(coef_sv)
                    coef_sh(ish) = sqrt(coef_psi(ish)*(coef_sv(ish)^2));
                end
                %     coef_sh = [(Sscale.sed_sv*sedv1),(Sscale.crust_sv*crustv1),(mantstructS1*mantlev1*mantlepsi1)];
            else
                coef_sv = coef_svn;
                coef_sh = coef_shn;
                
                coef_psi = (coef_sh.^2)./(coef_sv.^2);
            end
            
            coefi.v = coef_sv;
            coefi.psi = coef_psi;
            coefi.z = coef_z;
            
            TYPE = 'S';
            
            [sv,spsi,sz,coefS] = calc_spline_mod(spline,coefi,TYPE,ispsi);
            
            % Calculate SH from SV and psi
            sh = sqrt((sv.^2).*spsi);
            
            % Convert S-model to P-model
            clear coefi TYPE
            
            coef_pv = [(Pscale.sed_pv*sedv1),(Pscale.crust_pv*crustv1),(Pscale.mant_pv*mantlev1)];
            coef_ph = [(Pscale.sed_pv*sedv1),(Pscale.crust_pv*crustv1),(Pscale.mant_ph*mantlev1)];
            
            % Don't forget to introduce Vp/Vs
            coef_pv = coef_pv*vpvs;
            coef_ph = coef_ph*vpvs;
            
            coefi.v = coef_pv;
            coefi.h = coef_ph;
            coefi.z = coef_z;
            
            TYPE = 'P';
            
            [pv,ph,pz,coefP] = calc_spline_mod(spline,coefi,TYPE);
            
            
            %% Create Initial Model
            
            % Write model out to a MINEOS model card
            %             CARD = ['EARS_GS',num2str(ipara),'.card'];
            CARD = ['E',num2str(xi(ilat,ilon)),'_',num2str(yi(ilat,ilon)),'_',num2str(iter),'.card'];
            write_MINE_mod(bot,sedh1,crusth,sv,sh,sz,pv,ph,CARD);
            
            % Check the model if you want
            if ischeckmodel
                check_model(CARD);
            end
            
            % Run Mineos
            tic
            frun_mineos(CARD);
            toc
            return
            % Grab the estimated phase velocities
			clear TPHV TCPER SPHV SCPER            


			TYPE = 'T';
            [TCPER,TPHV,TUPER,TGRV,TMODE] = calc_fundCU(Cperiods,Uperiods,TYPE);
            TYPE = 'S';
            [SCPER,SPHV,SUPER,SGRV,SMODE] = calc_fundCU(Cperiods,Uperiods,TYPE);
            
            forward.tphv = TPHV;
            forward.tper = TCPER;
            forward.sphv = SPHV;
            forward.sper = SCPER;
            
            %% Run the inversion
%             mest = finverse_aniso(obs.Lphv,obs.Lstd,obs.Rphv,obs.Rstd,forward,kern);
            mest = finverse_aniso_linear(obs.Lphv,obs.Lstd,obs.Rphv,obs.Rstd,forward,frech_T,frech_S);
            
            %% update the model parameters
            % First SV
            ind_pos = find(mest.SV > maxpert);
            pertSV = mest.SV;
            pertSV(ind_pos) = maxpert;
            
            ind_neg = find(mest.SV < -1*maxpert);
            pertSV(ind_neg) = -1*maxpert;
            
            coef_svn = coef_sv+pertSV';
            
            % Next SH
            clear ind_pos
            ind_pos = find(mest.SH > maxpert);
            pertSH = mest.SH;
            pertSH(ind_pos) = maxpert;
            
            ind_neg = find(mest.SH < -1*maxpert);
            pertSH(ind_neg) = -1*maxpert;
            
            coef_shn = coef_sh+pertSH'; 
            % Perturb the p wave spline coefficients in accord with the s-wave
            coef_pvn = coef_pv.*(1+pertSV'./coef_sv);
            coef_phn = coef_ph.*(1+pertSH'./coef_sh);
            
            % save results 
            result(iter).forward.tphv = forward.tphv;
            result(iter).forward.tper = forward.tper;
            result(iter).forward.sper = forward.sper;
            result(iter).forward.sphv = forward.sphv;
            result(iter).card = CARD;
            result(iter).pertSV = pertSV;
            result(iter).pertSH = pertSH;
            result(iter).coef_sh = coef_sh;
            result(iter).coef_sv = coef_sv;
            result(iter).coef_svn = coef_svn;
            result(iter).coef_shn = coef_shn;

            
            figure(77)
            hold on
            plot(forward.tper,forward.tphv,'-k','color',CC(iter,:))
            plot(forward.sper,forward.sphv,'-k','color',CC(iter,:))
		            
            iter = iter+1;
            return
        end % max iter
        figure(77)
        hold on
        plot(Cperiods,obs.Lphv,'ok','linewidth',2)
        plot(Cperiods,obs.Rphv,'ok','linewidth',2)
%         result.obs.Lphv = obs.Lphv;
%         result.obs.Lstd = obs.Lstd;
%         result.obs.Rphv = obs.Rphv;
%         result.obs.Rstd = obs.Rstd;
        % plot the results of hte inversion
        return
        
    end % ilon
end % ilat
