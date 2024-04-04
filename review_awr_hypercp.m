% Input data from make_database_hypercp.m
% Add AVW and QWIP if not already included in Derived Products
% Display spectral (ir)radiance and reflectance with ancillary data
% dashboard

% Screen on QWIP,RelAz,Wind,SZA, and visual inspection
%
% D. Aurin NASA/GSFC March 2024

%% Setup
wipe
fontName = machine_prefs;

% Set to true unless building .env.all
validation = 0;

cruise = 'EXPORTSNA_NASA';
% cruise = 'EXPORTSNA_Boss';

clobber = 1;
plotQWIP = 0;
plotFlags = 0;
manualSelection = 0;

tRelAz = [88 137]; % M99, Z17, IOCCG
tSZA = [18 62]; % e.g. 20: Zhang 2017, depends on wind, e.g. 60:Brewin 2016
tWind = 10; %  6-7 m/s: IOCCG Draft Protocols, D'Alimonte pers. comm. 2019; 10 m/s: NASA SeaWiFS Protocols; 15 m/s: Zibordi 2009,
tQWIP = 0.2; % Dierssen et al. 2022
tQA = 0.2; % This is more experimental. Monitor it, but don't filter it.
tCloud = [20 80]; % Clear and fully overcast should be okay. 20% - 80% are likely to be the worst. This is experimental.
cloudIndexes = [0.05 0.3]; % From Ruddick et al. 2006 based on M99 models, where <0.05 is clear, >0.3 is fully overcast
if validation
    % Thresholds for validation
    tWind = 6.5; %  6-7 m/s: IOCCG Draft Protocols, D'Alimonte pers. comm. 2019; 10 m/s: NASA SeaWiFS Protocols; 15 m/s: Zibordi 2009,
    tQWIP = 0.17; % Experimental
end
%% Load and overview
if clobber
    load(sprintf('dat/%s.mat',cruise)) % SASData from make_database_hypercp.m
    dateTime = datetime([SASData.Ancillary.datenum],'ConvertFrom','datenum','TimeZone','UTC');
    nEns = size([SASData.Ancillary.datenum],2);

    try
        % Potential of wave mismatch depending on glint processing
        wave = vertcat(SASData.Reflectance(1).Rrs_HYPER_wavelength);
        Rrs = vertcat(SASData.Reflectance.Rrs_HYPER);
        Es = vertcat(SASData.Irradiance.ES_HYPER);
        Li = vertcat(SASData.Radiance.LI_HYPER);
        Lw = vertcat(SASData.Radiance.LW_HYPER);
    catch
        disp('Wavelength cannot be concatonated. Interpolating to shorter wave range')
        % In this case, interpolation of data to the shorter wavelength vector
        % will be required
        % There has to be a better way...
        for i=nEns:-1:1
            nBands(i) = length(SASData.Reflectance(i).Rrs_HYPER_wavelength);
        end
        iShort = find(nBands == min(nBands)); iShort = iShort(1);
        wave = SASData.Reflectance(iShort).Rrs_HYPER_wavelength;
        for i=nEns:-1:1
            Rrs(i,:) = interp1(SASData.Reflectance(i).Rrs_HYPER_wavelength, ...
                SASData.Reflectance(i).Rrs_HYPER, wave);
            Es(i,:) = interp1(SASData.Irradiance(i).ES_HYPER_wavelength, ...
                SASData.Irradiance(i).ES_HYPER, wave);
            Li(i,:) = interp1(SASData.Radiance(i).LI_HYPER_wavelength, ...
                SASData.Radiance(i).LI_HYPER, wave);
            Lw(i,:) = interp1(SASData.Radiance(i).LW_HYPER_wavelength, ...
                SASData.Radiance(i).LW_HYPER, wave);
        end
    end
    %% Apply the QWIP and AVW, regardless of whether they are provided
    [qwip,QCI,avw] = AVW_QWIP_2D_fun(Rrs,wave,'none','none');
    if plotQWIP
        minMaxPlot = [440 590];
        fh12 = QWIP_figure_fun(Rrs,wave,avw,cruise,minMaxPlot);
        exportgraphics(fh12(1),sprintf('plt/%s_QWIP.png',cruise))
        exportgraphics(fh12(2),sprintf('plt/%s_QWIP_Hist.png',cruise))
    end

    %% Populate threshold variables
    if sum(contains(fieldnames(SASData),'Derived_Products')) > 0
        if sum(contains(fieldnames(SASData.Derived_Products),'chlor_a')) > 0
            chl = [SASData.Derived_Products.chlor_a];
        else
            % Calculate chl
            chl=nan;
        end
    else
        % Calculate chl
        chl=nan;
    end
    if sum(contains(fieldnames(SASData),'Derived_Products')) > 0
        if sum(contains(fieldnames(SASData.Derived_Products),'QA_score')) > 0
            qa = [SASData.Derived_Products.QA_score];
        else
            % Calculate QA_score
            qa=nan;
        end
    else
        % Calculate QA_score
        qa=nan;
    end
    if sum(contains(fieldnames(SASData.Ancillary),'CLOUD')) > 0
        cloud = [SASData.Ancillary.CLOUD]; % Percent
    else
        cloud = nan(1,nEns);
    end
    % Fill in with cloud index where empty
    whrNaN = isnan(cloud)';
    [wv,iwv] = find_nearest(750,wave);
    li750 = Li(:,iwv);
    es750 = Es(:,iwv);
    cloud(whrNaN & (li750./es750 < cloudIndexes(1))) = 0; % Clear
    cloud(whrNaN & ( li750./es750 >= cloudIndexes(1)) & li750./es750 < cloudIndexes(2) ) = 50; % Partly cloudy
    cloud(whrNaN & (li750./es750 > cloudIndexes(2))) = 100; % Fully overcast (no good for validation, but okay for AWR)

    if sum(contains(fieldnames(SASData.Ancillary),'WINDSPEED')) > 0
        wind = [SASData.Ancillary.WINDSPEED];
    end
    if sum(contains(fieldnames(SASData.Ancillary),'SZA')) > 0
        sza = [SASData.Ancillary.SZA];
    end
    if sum(contains(fieldnames(SASData.Ancillary),'REL_AZ')) > 0
        relAz = [SASData.Ancillary.REL_AZ];
    end

    %% Set flags
    if validation
        flags.Cloud = cloud > tCloud(1) & cloud < tCloud(2); %Partly cloudy
    else
        flags.Cloud = false(1,length(cloud));  
    end
    flags.Wind = wind > tWind;
    flags.SZA = sza < tSZA(1) | sza > tSZA(2);
    flags.RelAz = abs(relAz)<tRelAz(1) | abs(relAz)>tRelAz(2);
    flags.QWIP = qwip > tQWIP;
    flags.QA = qa < tQA;

    save(sprintf('dat/%s_flags.mat',cruise),"dateTime","wave","Rrs","Es","Li","Lw",...
        "chl","avw","qwip","qa","cloud","wind","sza","relAz","flags")

    %% Apply flags
    flag = [flags.Cloud] | [flags.Wind] | [flags.SZA] | [flags.RelAz] | [flags.QWIP];

    %% Write CSV file for awr2env.py
    % In effect, we will have 2 files. One will have 0 or 1, the other 0 or 2.
    m(:,1) = year(dateTime);
    m(:,2) = month(dateTime);
    m(:,3) = day(dateTime);
    m(:,4) = hour(dateTime);
    m(:,5) = minute(dateTime);
    m(:,6) = second(dateTime);
    if validation
        csvOutFile = sprintf('dat/%s_flags.csv',cruise);
        m(:,7) = 2*int8(~flag); % 0, 1, or 2 for reject, seabass-only, validation
    else
        csvOutFile = sprintf('dat/%s_all_flags.csv',cruise);
        m(:,7) = int8(~flag); % 0, 1, or 2 for reject, seabass-only, validation
    end

    writematrix(m,csvOutFile)
else
    load(sprintf('dat/%s_flags.mat',cruise))
    nEns = size(Rrs,1);
end

%% Show spectra with filtered data
if plotFlags
    fh3 = figure;
    th1 = tiledlayout(2,2);
    ax1 = nexttile;
    plot(wave,Rrs,'k')
    hold on
    flagSpectra(ax1,wave,Rrs,flags,1)
    ylabel('R_{rs} [sr^{-1}]')

    ax2 = nexttile;
    plot(wave,Es,'k')
    hold on
    flagSpectra(ax2,wave,Es,flags,0)
    ylabel('E_s [\muW cm^{-2} nm^{-1}]')

    ax3 = nexttile;
    plot(wave,Li,'k')
    hold on
    flagSpectra(ax3,wave,Li,flags,0)
    ylabel('L_i [\muW cm^{-2} nm^{-1} sr^{-1}]')

    ax4 = nexttile;
    plot(wave,Lw,'k')
    hold on
    flagSpectra(ax4,wave,Lw,flags,0)
    ylabel('L_w [\muW cm^{-2} nm^{-1} sr^{-1}]')

    th1 = text(0.75,0.9,sprintf('Cloud: %d',sum(flags.Cloud)),'Units','normalized');
    th2 = text(0.75,0.85,sprintf('Wind: %d',sum(flags.Wind)),'Units','normalized');
    th3 = text(0.75,0.8,sprintf('SZA: %d',sum(flags.SZA)),'Units','normalized');
    th4 = text(0.75,0.75,sprintf('RelAz: %d',sum(flags.RelAz)),'Units','normalized');
    th5 = text(0.75,0.7,sprintf('QWIP: %d',sum(flags.QWIP)),'Units','normalized');
    th6 = text(0.75,0.65,sprintf('QA: %d',sum(flags.QA)),'Units','normalized','Color','r');

    gud = ~flags.RelAz & ~flags.Wind & ~flags.SZA & ~flags.QWIP & ~flags.Cloud; %& ~flags.QA
    th7 = text(0.75,0.6,sprintf('Remaining: %d of %d',sum(gud),nEns),'Units','normalized');

    set([ax1 ax2 ax3 ax4],'FontName',fontName,'FontSize',16, 'xgrid','on', 'ygrid','on')
    set(fh3,'position',[1926         381        1196         979])

    exportgraphics(fh3,sprintf('plt/%s_AllSpec.png',cruise))

    %% Figure timeline Chl, AVW, QWIP, Wei score, ..
    fh4 = figure;
    set(fh4,'Position',[200 200 850 950])
    ax1 = subplot(4,1,1);
    title(strrep(cruise,'_','-'))
    yyaxis(ax1,'left')
    ph1 = plot(dateTime,chl,'marker','.','markersize',18);
    ylabel('chlor_a')

    yyaxis(ax1,'right')
    ph2 = plot(dateTime,avw,'marker','.','markersize',18);
    ylabel('AVW [nm]')

    ax2 = subplot(4,1,2);
    yyaxis(ax2,'left')
    ph3 = plot(dateTime,qwip,'marker','.','markersize',18);
    ylabel('QWIP')

    yyaxis(ax2,'right')
    ph4 = plot(dateTime,qa,'marker','.','markersize',18);
    ylabel('QA_score')

    ax3 = subplot(4,1,3);
    yyaxis(ax3,'left')
    ph5 = plot(dateTime,cloud,'marker','.','markersize',18);
    ylabel('Cloud [%]')

    yyaxis(ax3,'right')
    ph6 = plot(dateTime,wind,'marker','.','markersize',18);
    ylabel('Wind [m/s]')

    ax4 = subplot(4,1,4);
    yyaxis(ax4,'left')
    ph7 = plot(dateTime,sza,'marker','.','markersize',18);
    ylabel('SZA')

    yyaxis(ax4,'right')
    ph8 = plot(dateTime,relAz,'marker','.','markersize',18);
    ylabel('Relative Azimuth')

    set([ax1 ax2 ax3 ax4],'fontname',fontName, 'fontsize', 14, 'xgrid', 'on')
    if validation
        exportgraphics(fh4,sprintf('plt/%s_timeline.png',cruise))
    else
        exportgraphics(fh4,sprintf('plt/%s_all_timeline.png',cruise))
    end


    %% Manual screening of spectra
    if manualSelection
        input('Continue now? (enter)');
        close all
        fh5 = figure;
        set(fh5,'position',[1926         381        1196         979])
        plot(wave,Rrs,'k')
        hold on
        flagSpectra(ax1,wave,Rrs,flags,0)
        ylabel('R_{rs} [sr^{-1}]')

        disp('Manual Screening')
        disp('Zoom to spectrum of interest and hit Continue')
        disp('Left mouse to flag spectrum. Continue to save and move on. Try again to ignore last.')
        disp('Middle mouse to exit and save.')
        flag = zeros(1,nEns);
        for tries=1:10
            h1 = uicontrol(fh5,'Style', 'pushbutton', 'String', 'Continue',...
                'Position', [5 100 50 25], 'Callback', 'uiresume');
            uicontrol(h1)
            uiwait(fh5)

            [x,y] = ginput(1);

            plot(x,y,'k*')
            h2 = uicontrol('Style', 'pushbutton', 'String', 'Continue',...
                'Position', [5 100 50 25], 'Callback', 'butt=1; uiresume');
            h3 = uicontrol('Style', 'pushbutton', 'String', 'Try Again',...
                'Position', [5 50 50 25], 'Callback','butt=0; uiresume');
            h4 = uicontrol('Style', 'pushbutton', 'String', 'Exit',...
                'Position', [5 5 50 25], 'Callback','butt=3; uiresume');
            uicontrol(h2)
            uicontrol(h3)
            uicontrol(h4)
            uiwait
            if butt == 0
                continue
            elseif butt == 1
                [wv,windex] = find_nearest(x,wave);
                RrsX = Rrs(:,windex);
                [rrs,Rindex] = find_nearest(y,RrsX);
                plot(wave,Rrs(Rindex,:),'k','LineWidth',3)
                flag(Rindex) = 1;
                fprintf('Index: %d\n',Rindex)
                disp(flag(Rindex))
                % disp([x,y])
                % continue
            elseif butt == 3
                break
            end
        end
    end
    flags.Manual = flag;
    %% Apply flags
    flag = [flags.Cloud] | [flags.Wind] | [flags.SZA] | [flags.RelAz] | [flags.QWIP] |...
        [flags.Manual];

    %% Write CSV file for awr2env.py
    % In effect, we will have 2 files. One will have 0 or 1, the other 0 or 2.
    if validation
        csvOutFile = sprintf('dat/%s_flags.csv',cruise);
        m(:,7) = 2*int8(~flag); % 0, 1, or 2 for reject, seabass-all, validation
    else
        csvOutFile = sprintf('dat/%s_all_flags.csv',cruise);
        m(:,7) = int8(~flag); % 0, 1, or 2 for reject, seabass-all, validation
    end

    m(:,1) = year(dateTime);
    m(:,2) = month(dateTime);
    m(:,3) = day(dateTime);
    m(:,4) = hour(dateTime);
    m(:,5) = minute(dateTime);
    m(:,6) = second(dateTime);

    writematrix(m,csvOutFile)
end

save(sprintf('dat/%s_flags.mat',cruise),"dateTime","wave","Rrs","Es","Li","Lw",...
    "chl","avw","qwip","qa","cloud","wind","sza","relAz","flags")


%%
function flagSpectra(ax,wave,Var,flags,leg)

if sum(flags.Cloud)>0
    ph1 = plot(ax,wave,Var(flags.Cloud,:),'y','linewidth',1);
else
    ph1 = plot(ax,wave,Var*nan,'y','linewidth',1);
end
if sum(flags.Wind)>0
    ph2 = plot(ax,wave,Var(flags.Wind,:),'r','linewidth',2);
else
    ph2 = plot(ax,wave,Var*nan,'r','linewidth',2);
end
if sum(flags.SZA)>0
    ph3 = plot(ax,wave,Var(flags.SZA,:),'g','linewidth',2);
else
    ph3 = plot(ax,wave,Var*nan,'g','linewidth',2);
end
if sum(flags.RelAz)>0
    ph4 = plot(ax,wave,Var(flags.RelAz,:),'b','linewidth',2);
else
    ph4 = plot(ax,wave,Var*nan,'b','linewidth',2);
end
if sum(flags.QWIP)>0
    ph5 = plot(ax,wave,Var(flags.QWIP,:),'m','linewidth',2);
else
    ph5 = plot(ax,wave,Var*nan,'m','linewidth',2);
end
if sum(flags.QA)>0
    ph6 = plot(ax,wave,Var(flags.QA,:),'b','linewidth',2,'linestyle','--');
else
    ph6 = plot(ax,wave,Var*nan,'b','linewidth',2,'linestyle','--');
end


% These won't concatonate
% set([ph1 ph2 ph3 ph4 ph5 ph6],'linewidth',2)
if leg
    legend([ph1(1) ph2(1)  ph3(1) ph4(1)  ph5(1)  ph6(1)],...
        fieldnames(flags))
end

end