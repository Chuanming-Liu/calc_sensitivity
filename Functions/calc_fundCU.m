% This is hardwired when in the spheroidal case because the card files are
% named with set PERIODS (e.g. 0_30, 31_70, 71_120)
% Calls the read_forward_model.m function to read the ascii cards output
% from MINEOS;
% Returns the period and phase velocity closest to the phase velocities of
% interest
%
% NJA, 2014
%
% Dont be scared by the use of the card name from the parameter file ... we
% are overwriting the ascii files from each run every time thus we only
% need to read the same card file
%

function [CPER,PHV,UPER,GRV,MODE] = calc_fundCU(CPERIODS,UPERIODS,TYPE)

setup_parameters;

CARD = param.CARDID;
isfigure = 0;

% % Turn on if want to run as a standalone script 
% UPERIODS = param.Uperiods;
% CPERIODS = param.Cperiods;
% TYPE = 'S';


if strcmp(TYPE,'T') == 1
    
    Tmode = read_forward_model(CARD,TYPE);
    
    fundind = find(Tmode.n == 0);
    
    fund_grv = Tmode.grv(fundind);
    fund_t = Tmode.t(fundind);
    fund_phv = Tmode.phv(fundind);
    
    % Phase Velocity
    for icp = 1:length(CPERIODS)
        match = abs(fund_t - CPERIODS(icp));
        tind = find(match == min(match));
                if length(tind) > 1
            for itind = 1:length(tind)
                disp(['PERIOD : ',num2str(fund_t(tind(itind))),' PHASEV : ',num2str(fund_phv(tind(itind)))]);
            end
            tind = tind(end);
        end
        CPER(icp) = fund_t(tind);
        PHV(icp) = fund_phv(tind);
    end
	% disp(CPER)
	clear tind match
    
    % Group Velocity
    for iup = 1:length(UPERIODS)
        clear tind match
        match = abs(fund_t - UPERIODS(iup));
        tind = find(match == min(match));
        if length(tind) > 1
            for itind = 1:length(tind)
                disp(['PERIOD : ',num2str(fund_t(tind(itind))),' PHASEV : ',num2str(fund_phv(tind(itind)))]);
            end
            tind = tind(end);
        end
        UPER(iup) = fund_t(tind);
        GRV(iup) = fund_grv(tind);
    end
    
    if isfigure
        figure(3)
        clf
        plot(fund_t,fund_grv,'-r','linewidth',2)
        hold on
        plot(fund_t,fund_phv,'-b','linewidth',2)
        plot(CPER,PHV,'ok','linewidth',2);
        plot(UPER,GRV,'ok','linewidth',2);
%         xlim([min(CPERIODS)-5 max(CPERIODS)+5])
        ylim([2.5 5])
        xlim([min(UPERIODS)-5 max(CPERIODS)+5])
        set(gca,'fontsize',16);
        title('Fundamental Mode');
    end
    
    % return the entire Tmode structure in case we want to look at things
    % later on
    MODE = Tmode;
    
elseif strcmp(TYPE,'S') == 1
    
    TYPE_long = 'S0_50';
    TYPE_short = 'S50_150';
    Smode1 = read_forward_model(CARD,TYPE_long);
    Smode2 = read_forward_model(CARD,TYPE_short);
    
    % Combine the two mode files together
    Smode.n = [Smode1.n',Smode2.n'];
    Smode.l = [Smode1.l',Smode2.l'];
    Smode.wrad = [Smode1.wrad',Smode2.wrad'];
    Smode.w = [Smode1.w',Smode2.w'];
    Smode.t = [Smode1.t',Smode2.t'];
    Smode.grv = [Smode1.grv',Smode2.grv'];
    Smode.phv = [Smode1.phv,Smode2.phv];
    Smode.rayq = [Smode1.rayq',Smode2.rayq'];
    Smode.q = [Smode1.q',Smode2.q'];
    
    % Pull out fundamental mode
    fundind = find(Smode.n == 0);
    
    fund_grv = Smode.grv(fundind);
    fund_t = Smode.t(fundind);
    fund_phv = Smode.phv(fundind);
    
    % Phase Velocity
    for icp = 1:length(CPERIODS)
        match = abs(fund_t - CPERIODS(icp));
        tind = find(match == min(match));
        if length(tind) > 1
            for itind = 1:length(tind)
                disp(['PERIOD : ',num2str(fund_t(tind(itind))),' PHASEV : ',num2str(fund_phv(tind(itind)))]);
            end
            tind = tind(end);
        end
        CPER(icp) = fund_t(tind);
        PHV(icp) = fund_phv(tind);
    end
    
    % Group Velocity
    for iup = 1:length(UPERIODS)
        match = abs(fund_t - UPERIODS(iup));
        tind = find(match == min(match));
                if length(tind) > 1
            for itind = 1:length(tind)
                disp(['PERIOD : ',num2str(fund_t(tind(itind))),' PHASEV : ',num2str(fund_phv(tind(itind)))]);
            end
            tind = tind(end);
        end
        UPER(iup) = fund_t(tind);
        GRV(iup) = fund_grv(tind);
    end
    
    if isfigure
        figure(1)
        clf
        plot(fund_t,fund_grv,'-b','linewidth',2)
        hold on
        plot(fund_t,fund_phv,'--r','linewidth',2)
        plot(CPER,PHV,'.k','markersize',30);
        plot(UPER,GRV,'or','linewidth',2);
        xlim([min(UPERIODS)-5 max(CPERIODS)+5])
        ylim([2.5 5])
        set(gca,'fontsize',16);
        title('Fundamental Mode');
    end
    
    % Return entire S mode structure in case we want to look at it later
    MODE = Smode;
else
    disp('Type not reqcognized!');
end
