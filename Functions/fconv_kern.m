% Perturb the spline model and calculation the perturbation to the
% individual velocity layers
% REMEMBER --- This is all interms of m/s not km/s!!!
% NJA, 4-29-2014
%
% function [DA] = perturb_mod(sv,sh,sz,coef,spline)
%
% 9/27/2014 - modified to output SV and SH sensitivities for both wave
% types
%
% clear

function [kern] = fconv_kern(periods,CARD,FRECH_T,FRECH_S,spline,coef_s)

isfigure = 0;

% bot = param.bot;
% crusth = param.mohoh;
% sedh = param.sedh;
% 
% periods = param.periods;
% CARD = param.CARDID;

% Load information about spline model
% load([CARD,'_iter0.mat']);

% spline = forward.spline;
% coef_s = forward.coef_s;

coef_sv = coef_s.sv;
coef_sh = coef_s.sh;
coef_z = coef_s.z;

% Load information about phase velocity frechet kernels
% load([CARD,'_fcv.mat']);

%% Initialize vectors for the splines
z_bs = [spline(1).z,spline(2).z,spline(3).z];

% Original velocities
sv_os = ones(length(z_bs),1);
sv_ns = ones(length(z_bs),1);
sh_os = ones(length(z_bs),1);
sh_ns = ones(length(z_bs),1);
%% Loop through splines

% [dum nlayer] = size(spline);

% Loop through layer - sed, crust, mant
iz = 1; % keep count of depth
ic = 1; % keep count of coefficients

dCb = 10; % perturbation to spline coefficients (Cb)


for itype = 1:2 % T and S
    clear FRECH dAdC
    figure(85)
    clf
    figure(83)
    clf
    figure(82)
    clf
    
    if itype == 1; FRECH = FRECH_T; % Toroidal
    elseif itype == 2; FRECH = FRECH_S; % Spheroidal
    end
    
    clear depth sv_coef sh_coef N
    
    %
    % Perturb each coefficient one by one
    DC = jet(length(coef_sv));
    for idc = 1:length(coef_sv)
        
        sv_coefn = coef_sv;
        sv_coefn(idc) = sv_coefn(idc)+dCb;
        
        sh_coefn = coef_sh;
        sh_coefn(idc) = sh_coefn(idc)+dCb;
        
        % Create new vector of coefficient to pass to function
        coef_new.v = sv_coefn;
        coef_new.h = sh_coefn;
        coef_new.z = coef_z;
        
        coef_old.v = coef_s.sv;
        coef_old.h = coef_s.sh;
        coef_old.z = coef_s.z;
        
        % calculate velocity model given new coefficient
        
        [sv_ns,sh_ns,sz_ns] = calc_spline_mod(spline,coef_new,'S');
        
        [sv_os,sh_os,sz_os] = calc_spline_mod(spline,coef_old,'S');
        
        %         if isfigure
        %             figure(81)
        %             clf
        %             hold on
        %             plot(sv_os,sz_os,'-k','linewidth',2)
        %             plot(sv_ns,sz_ns,':r','linewidth',2)
        %             plot(coef_new.sv,coef_new.z,'or','linewidth',2);
        %             plot(coef_old.sv,coef_old.z,'ok','linewidth',2);
        %             xlim([1000 4500])
        %             ylim([0 200])
        %             set(gca,'ydir','reverse')
        %             pause
        %         end
        
        %% Now calculate dB/dC
        % Determine which layers were affected
        clear dsv dsh;
        dsv = sv_ns-sv_os;
        dsh = sh_ns-sh_os;
        
        dBsvdC = zeros(size(dsv));
        dBshdC = zeros(size(dsh));
        dAdSV = zeros(size(dsv));
        dAdSH = zeros(size(dsv));
        dAdSVdBdC = zeros(size(dsv));
        dAdSHdBdC = zeros(size(dsv));
        
        dC = dCb;
        dBshdC = dsh/dC;
        dBsvdC = dsv/dC;
        
        % Now load in sensitivity information from frechet_cv file
        % For now we only care about Vsv and Vsh
        CC = jet(length(periods));
        
        % initialize count variables to find where the problem is occuring
        count1 = 0;
        count2 = 0;
        for ip = 1:length(periods)
            fsv = FRECH(ip).vsv;
            fsh = FRECH(ip).vsh;
            fz = 6371-FRECH(ip).rad./1000;
            
            % Now run through all depths within the frechet file and find
            % the correct dBdC value
            
            for iz = 1:length(sz_ns)
                vind = find(fz == sz_ns(iz));
                
                % Check to see if more than one depth matches
                % -- This is actually not an issue as the locations where
                % this happens is at the surface (z = 0) and at the
                % sediments and moho where we expect duplications to occur
                % in the frechet file.
                if isempty(vind) == 1
                    %                     disp('No depth found!')
                    %                     count1 = count1+1;
                    %                     zerocount(count1) = sz_ns(iz);
                    continue
                elseif length(vind) >1
                    %                     disp('More than one depth matches -- Likely a discontinuity');
                    %                     count2 = count2 +1;
                    %                     multcount(count2,1) = fz(vind(1));
                    %                     multcount(count2,2) = sz_ns(iz);
                    vind = vind(end);
                end;
                
                % fprintf('FZ : %f SZ : %f \n',[fz(vind) sz_ns(iz)]);
                
                dAdSV(iz) = fsv(vind);
                dAdSH(iz) = fsh(vind);
                
                dAdSVdBdC(iz) = dAdSV(iz)*dBsvdC(iz);
                dAdSHdBdC(iz) = dAdSH(iz)*dBshdC(iz);
                
            end
            
            dASHdC(idc,ip) = sum(dAdSHdBdC); % idc is the index of the coefficient, ip for the period of interest
            dASVdC(idc,ip) = sum(dAdSVdBdC);
            
            
            
            if isfigure;
                % Compare kernels read from the file and the same kernels
                % but with our indexing
                figure(83)
                subplot(1,2,1)
                hold on
                %                 plot(fsv,fz,'-k','linewidth',2,'color',CC(ip,:))
                plot(fsv,fz,'-k','linewidth',2,'color','b')
                plot(dAdSV,sz_ns,'-k','linewidth',1);
                ylim([0 200])
                set(gca,'ydir','reverse','fontsize',14)
                title('SV')
                subplot(1,2,2)
                hold on
                %                 plot(fsh,fz,'-k','linewidth',2,'color',CC(ip,:));
                plot(fsh,fz,'-k','linewidth',2,'color','r');
                
                plot(dAdSH,sz_ns,'-k','linewidth',1);
                ylim([0 200])
                set(gca,'ydir','reverse','fontsize',14)
                title('SH')
            end;
            
        end % end loop over periods
        
        figure(85)
        subplot(1,2,1)
        hold on
        plot(periods,dASVdC(idc,:),'-k','linewidth',2,'color',DC(idc,:))
        set(gca,'fontsize',14)
        xlim([20 110]);
        xlabel('PERIOD (S)')
        ylabel('dA/dC')
        title('SV')
        subplot(1,2,2)
        hold on
        plot(periods,dASHdC(idc,:),'-k','linewidth',2,'color',DC(idc,:))
        xlim([20 110]);
        set(gca,'fontsize',14)
        xlabel('PERIOD (S)')
        ylabel('dA/dC')
        title('SH')
        
        if isfigure
            figure(82)
            subplot(1,2,1)
            hold on
            plot(dsv,sz_ns,'-k','color',DC(idc,:),'linewidth',2)
            set(gca,'ydir','reverse')
            title('SV')
            subplot(1,2,2)
            hold on
            plot(dsh,sz_ns,'-k','color',DC(idc,:),'linewidth',2)
            set(gca,'ydir','reverse')
            title('SH')
            
        end
    end % end loop over coefficients
    
    if itype == 1 % toroidal
        kern.TSV = dASVdC;
        kern.TSH = dASHdC;
    elseif itype == 2 % spheroidal
        kern.SSV = dASVdC;
        kern.SSH = dASHdC;
    end
    
    
end % end loop over T or S

% savefile = [CARD,'_kern.mat'];

% save(savefile,'kern');
