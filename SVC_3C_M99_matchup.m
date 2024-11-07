% Load database from cruises with both 3C and Mobley rho corrections and
% compare the two
%
% Input:
%   dBase structure from make_awr_seabass.m
%
% Output:
%   Plots of matched stations 3C vs. M99
%
% D. Aurin NASA/GSFC November 2024

wipe
fontName = machine_prefs;
% cruise = 'PVST_PRINGLS_PRINGLS_20240417';
cruise = 'PVST_PRINGLS_PRINGLS_20240513';

load(sprintf('dat/%s.mat',cruise)) % dBase from make_awr_seabass.m

rho = extractfield(dBase,'rho_correction');
station = extractfield(dBase,'station');

C3ind = strcmp(rho,'3C');
M99ind = strcmp(rho,'M99');

C3 = dBase(C3ind);
M99 = dBase(M99ind);

plotSpec = 1;
plotScat = 1;
plotBA = 1;
plotMet = 1;

%% Spectral
if plotSpec
    cmap = jet(length(C3));
    fh1 = figure;
    for i=1:length(C3)
        if i==1
            ph1 = plot(C3(i).wavelength,C3(i).rrs,'Color',cmap(i,:));
            hold on
            ph2 = plot(M99(i).wavelength,M99(i).rrs,'Color',cmap(i,:),'Marker','.','LineStyle','none');
        else
            plot(C3(i).wavelength,C3(i).rrs,'Color',cmap(i,:))
            plot(M99(i).wavelength,M99(i).rrs,'Color',cmap(i,:),'Marker','.','LineStyle','none')
        end
    end
    grid on
    xlabel('wavelength [nm]','FontName',fontName,'Fontsize',14)
    ylabel('R_{rs} [sr^{-1}]','FontName',fontName,'Fontsize',14)
    legend([ph1 ph2],'3C','M99','Fontsize',18)
    set(gca,'position',[0.1300    0.1100    0.7750    0.8150],'FontName',fontName,'Fontsize',14)
    exportgraphics(fh1,sprintf('plt/3C_M99_spectral_%s.png',cruise))
end

%% Scatter
bands = [380 410 550 620 700 800];
wv = C3(1).wavelength;
% wavelengths appear same for all stations, both rho corrections    
rrs3C = vertcat(C3.rrs);
rrsM99 = vertcat(M99.rrs);
if plotScat        
    fh2 = figure;
    for i=1:length(bands)
        [~,ind] = find_nearest(bands(i),wv);
        subplot(3,2,i)
        plot(rrs3C(:,ind),rrsM99(:,ind),'.','MarkerSize',14,'Color','r')
        p11
        grid on
        title(bands(i))
        ylabel('M99')
        xlabel('3C')
        set(gca,'FontName',fontName,'Fontsize',14)
    end
    set(gcf,'Position',[1841         205         826         715])
    exportgraphics(fh2,sprintf('plt/3C_M99_scatter_%s.png',cruise))
end
close all

%% Bland-Altman
if plotBA
    for i=1:length(bands)
        [~,ind] = find_nearest(bands(i),wv);
        % subplot(3,2,i)
        names = {sprintf('3C %d',bands(i)),sprintf('M99 %d',bands(i))};
        [rpc fh3] = BlandAltman(rrs3C(:,ind),rrsM99(:,ind),names);%,'.','MarkerSize',14,'Color','r')
        exportgraphics(fh3,sprintf('plt/3C_M99_BA_%d_%s.png',bands(i),cruise))
        close
    end
end

%% Wind/Cloud
if plotMet
    wind = [C3.wind];
    cloud = [C3.cloud];
    fh4 = figure;
    for i=1:length(bands)
        [~,ind] = find_nearest(bands(i),wv);
        subplot(3,2,i)
        yyaxis left
        plot(rrsM99(:,ind)-rrs3C(:,ind),wind,'.','MarkerSize',14);%,'Color','r')
        ylabel('Wind')
        yyaxis right
        plot(rrsM99(:,ind)-rrs3C(:,ind),cloud,'.','MarkerSize',14);
        ylabel('Cloud')

        grid on
        title(bands(i))
        xlabel('M99 - 3C')

        set(gca,'FontName',fontName,'Fontsize',14)
    end
    set(gcf,'Position',[1841         205         826         715])
    exportgraphics(fh4,sprintf('plt/3C_M99_wind_cloud_%s.png',cruise))
end




