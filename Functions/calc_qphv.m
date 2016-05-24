% Function to extract q corrected phase velocity from mineos
%
% NJA, 1/31/2014

function [qcor,ucor] = calc_qphv(periods)

setup_parameters;
datapath = param.DATAPATH;
SID = param.SID;
TID = param.TID;
CARDID = param.CARDID;

isfigure = 0;
for it =1:2
    if it == 1
        % spheroidal first
        TYPE = 'S'
        QASC = [datapath,CARDID,'.',SID,'.q'];
        colr = 'g';
    elseif it == 2
        % toroidal next
        TYPE = 'T';
        QASC = [datapath,CARDID,'.',TID,'.q'];
        colr = 'r';
    end;
    
    % Paramters for reading the file
    hline = 12;
    
    % Read card information
    fid = fopen(QASC,'r');
    
    A = textscan(fid,'%d %d %f %f %f %f %f %f %f %f','headerlines',hline);
    
    qcard.name = QASC;
    qcard.n = A{1};
    qcard.l = A{2};
    qcard.wrad = A{3};
    qcard.qq = A{4};
    qcard.phi = A{5};
    qcard.phv = A{6};
    qcard.grv = A{7};
    qcard.qphv = A{8};
    qcard.Tq = A{9};
    qcard.T = A{10};
    
    subplot(1,2,it)
    hold on
    for ip = 1:length(periods)
        qind = find(abs(qcard.Tq - periods(ip)) == min(abs(qcard.Tq - periods(ip))));
        oind = find(abs(qcard.T - periods(ip)) == min(abs(qcard.T - periods(ip))));
        if isempty(oind)
            disp(['No original period found for ',num2str(periods(ip))])
            pause
        elseif isempty(qind)
            disp(['No q corrected period found for ',num2str(periods(ip))])
        end
        
        if isfigure
            if it == 1 && ip == 1
                figure(37)
                clf
            end
            plot(qcard.T(oind),qcard.phv(oind),'ok','linewidth',2,'markerfacecolor','k','markersize',10);
            plot(qcard.Tq(qind),qcard.qphv(qind),'ok','linewidth',2,'markerfacecolor',colr,'markersize',10)
            xlim([5 100])
        end
        if TYPE == 'T'
            qcor.tper(ip) = qcard.Tq(qind);
            qcor.tphv(ip) = qcard.qphv(qind);
            ucor.tper(ip) = qcard.T(oind);
            ucor.tphv(ip) = qcard.phv(oind);
        elseif TYPE == 'S'
            qcor.sper(ip) = qcard.Tq(qind);
            qcor.sphv(ip) = qcard.qphv(qind);
            ucor.sper(ip) = qcard.T(oind);
            ucor.sphv(ip) = qcard.phv(oind);
        end
    end
end

