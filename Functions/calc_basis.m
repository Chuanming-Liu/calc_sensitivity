
% Create cubic splines for model setup
% Thicknesses of layers are not allowed to change
% This script can be easily converted to a function that takes in the layer
% depths and gives back a structure of all the splines
% NJA, 2014
%
% Modified 9/2014 so that mantle knots defined by the bottom depth

function [spline] = calc_basis(bot,sedh,crusth);

isfigure = 0;
issaved = 0;

% % Turn on if running not as a function
% clear
% bot = 200;
% sedh = 5;
% crusth = 40;

% Calculate the Location of Knots in the Mantle
diff = bot-crusth;

mantk1 = crusth+(diff*.25);
mantk2 = crusth+(diff*.6);

% First the sediments -- linear
% z_sed = 0:0.5:sedh;
z_sed = 0:0.5:sedh;
t_sed = [0 sedh sedh];
bs_sed = bspline_basis(0,1,t_sed,z_sed);

% Now the crust -- 3 cubic
z_crust = sedh:crusth;

t_crust = [sedh sedh sedh crusth crusth crusth];
bs_crust1 = bspline_basis(0,3,t_crust,z_crust);
bs_crust2 = bspline_basis(1,3,t_crust,z_crust);
bs_crust3 = bspline_basis(2,3,t_crust,z_crust);


% Next the mantle -- 5 cubic
z_mant = crusth:bot;
t_mant = [crusth crusth crusth mantk1 mantk2 bot bot bot];
bs_mant1 = bspline_basis(0,3,t_mant,z_mant);
bs_mant2 = bspline_basis(1,3,t_mant,z_mant);
bs_mant3 = bspline_basis(2,3,t_mant,z_mant);
bs_mant4 = bspline_basis(3,3,t_mant,z_mant);
bs_mant5 = bspline_basis(4,3,t_mant,z_mant);

if isfigure
    figure(71)
    clf
    hold on
    box on
    set(gca,'LineWidth',2)
    
    plot(bs_sed,z_sed,'-g','linewidth',2);
    
    plot(bs_crust1,z_crust,'-r','linewidth',2);
    plot(bs_crust2,z_crust,'-r','linewidth',2);
    plot(bs_crust3,z_crust,'-r','linewidth',2);
%     plot(bs_crust4,z_crust,'-r','linewidth',2);
    
    plot(bs_mant1,z_mant,'-b','linewidth',2);
    plot(bs_mant2,z_mant,'-b','linewidth',2);
    plot(bs_mant3,z_mant,'-b','linewidth',2);
    plot(bs_mant4,z_mant,'-b','linewidth',2);
    plot(bs_mant5,z_mant,'-b','linewidth',2);
    
    ylim([0 bot])
    ylabel('Depth (km)')
    set(gca,'ydir','reverse')
    
    xdum = 0:.05:1;
    figure(72)
    clf
    hold on
    box on
    set(gca,'LineWidth',2)
        plot(bs_sed,z_sed,'-g','linewidth',2);
    
    plot(bs_crust1,z_crust,'-r','linewidth',2);
    plot(bs_crust2,z_crust,'-r','linewidth',2);
    plot(bs_crust3,z_crust,'-r','linewidth',2);
%     plot(bs_crust4,z_crust,'-r','linewidth',2);
    
    plot(bs_mant1,z_mant,'-b','linewidth',2);
    plot(bs_mant2,z_mant,'-b','linewidth',2);
    plot(bs_mant3,z_mant,'-b','linewidth',2);
    plot(bs_mant4,z_mant,'-b','linewidth',2);
    plot(bs_mant5,z_mant,'-b','linewidth',2);
    
    plot(xdum,ones(size(xdum))*mantk1,'.k','markersize',10);
    plot(xdum,ones(size(xdum))*mantk2,'.k','markersize',10);
        
    ylim([0 bot])
    ylabel('Depth (km)')
    set(gca,'ydir','reverse')
    
end

% Save everuthing to a mat file for safe keeping
spline(1).z = z_sed;
spline(2).z = z_crust;
spline(3).z = z_mant;

spline(1).t = t_sed;
spline(2).t = t_crust;
spline(3).t = t_mant;

spline(1).b(1,:) = bs_sed;

spline(2).b(1,:) = bs_crust1;
spline(2).b(2,:) = bs_crust2;
spline(2).b(3,:) = bs_crust3;

spline(3).b(1,:) = bs_mant1;
spline(3).b(2,:) = bs_mant2;
spline(3).b(3,:) = bs_mant3;
spline(3).b(4,:) = bs_mant4;
spline(3).b(5,:) = bs_mant5;


if issaved == 1;
    savefile = 'spline.mat';
    save(savefile,'spline');
end