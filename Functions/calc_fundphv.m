% Function to pick out estimated phase velocities from calculated
% dispersion curves.
% This is hardwired when in the spheroidal case because the card files are
% named with set PERIODS (e.g. 0_30, 31_70, 71_120)
% Calls the read_forward_model.m function to read the ascii cards output
% from MINEOS;
% Returns the period and phase velocity closest to the phase velocities of
% interest
%
% NJA, 2014

function [PER,PHV] = calc_fundphv(PERIODS,TYPE)

setup_parameters;


% PERIODS = param.PERIODS;
% TYPE = 'T';

CARD = param.CARDID;
isfigure = 0;

if strcmp(TYPE,'T') == 1
    Tmode = read_forward_model(CARD,TYPE);
    
    fundind = find(Tmode.n == 0);

    fund_grv = Tmode.grv(fundind);
    fund_t = Tmode.t(fundind);
    fund_phv = Tmode.phv(fundind);
    
            for ip = 1:length(PERIODS)
            match = abs(fund_t - PERIODS(ip));
            tind = find(match == min(match));
            
            PER(ip) = fund_t(tind);
            PHV(ip) = fund_phv(tind);
        end
    
   if isfigure
        figure(3)
        clf
        plot(fund_t,fund_grv,'-r','linewidth',2)
        hold on
        plot(fund_t,fund_phv,'--r','linewidth',2)
        plot(PER,PHV,'ok','linewidth',2);
        xlim([min(PERIODS)-5 max(PERIODS)+5])
        ylim([2.5 5])
        set(gca,'fontsize',16);
        title('Fundamental Mode');
    end
    
    
    
elseif strcmp(TYPE,'S') == 1
             
    Smode = read_forward_model(CARD,TYPE);
    
    fundind = find(Smode.n == 0);

    fund_grv = Smode.grv(fundind);
    fund_t = Smode.t(fundind);
    fund_phv = Smode.phv(fundind);
         
        for ip = 1:length(PERIODS)
            match = abs(fund_t - PERIODS(ip));
            tind = find(match == min(match));
            
            PER(ip) = fund_t(tind);
            PHV(ip) = fund_phv(tind);
        end

        if isfigure
        figure(1)
        clf
        plot(fund_t,fund_grv,'-b','linewidth',2)
        hold on
        plot(fund_t,fund_phv,'--b','linewidth',2)
        plot(PER,PHV,'ok','linewidth',2);
        xlim([min(PERIODS)-5 max(PERIODS)+5])
        ylim([2.5 5])
        set(gca,'fontsize',16);
        title('Fundamental Mode');
    end
else
    disp('Type not reqcognized!');
end





