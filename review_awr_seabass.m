% Input data from make_database_seabass.m
% Add AVW and QWIP if not already included in Derived Products
% Display spectral (ir)radiance and reflectance with ancillary data
% dashboard

% Screen on QWIP,RelAz,Wind,SZA, and visual inspection
%
% D. Aurin NASA/GSFC June 2024

%% Setup
wipe
fontName = machine_prefs;

% Set to true unless building .env.all
validation = 0;

cruise = 'UMCES_Missouri_Reservoirs';
SBA = 1;

clobber = 1;
plotQWIP = 0;
plotFlags = 1;
manualSelection = 1;

negRrs = [380 680]; % Spectral range of negatives to eliminate from all sets
tRelAz = [87 138]; % M99, Z17, IOCCG
tSZA = [18 62]; % e.g. 20: Zhang 2017, depends on wind, e.g. 60:Brewin 2016
tWind = 10; %  6-7 m/s: IOCCG Draft Protocols, D'Alimonte pers. comm. 2019; 10 m/s: NASA SeaWiFS Protocols; 15 m/s: Zibordi 2009,
tQWIP = 0.2; % Dierssen et al. 2022
tQA = 0.2; % This is more experimental. Monitor it, but don't filter it.
tCloud = [20 80]; % Clear and fully overcast should be okay. 20% - 80% are likely to be the worst. This is experimental.
cloudIndexes = [0.05 0.3]; % From Ruddick et al. 2006 based on M99 models, where <0.05 is clear, >0.3 is fully overcast
if validation
    % Thresholds for validation
    tRelAz = [89 136]; % M99, Z17, IOCCG
    tWind = 7; %  6-7 m/s: IOCCG Draft Protocols, D'Alimonte pers. comm. 2019; 10 m/s: NASA SeaWiFS Protocols; 15 m/s: Zibordi 2009,
    tQWIP = 0.17; % Experimental
end
%% Load and overview
if clobber
    load(sprintf('dat/%s.mat',cruise)) % dBase from make_awr_seabass.m
    dateTime = [dBase.datetime];
    nEns = length(dBase);

    try
        % Potential of wave mismatch depending on glint processing
        wave = dBase(1).wavelength;
        Rrs = vertcat(dBase.rrs);
        Es = vertcat(dBase.es);
        Lw = vertcat(dBase.lw);
        if ~SBA
            % Have not encountered this yet
            Li = vertcat(dBase.li);            
        end
    catch
        disp('Wavelength cannot be concatonated. Interpolating to shorter wave range')
        % Cross this bridge when/if it arises
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
    if sum(contains(fieldnames(dBase),'chlor_a')) > 0
        chl = [dBase.chlor_a];
    else
        % Calculate chl
        chl=nan(1,nEns);
    end
    if sum(contains(fieldnames(dBase),'QA_score')) > 0
        qa = [dBase.QA_score];
    else
        % Calculate QA_score
        qa=nan(1,nEns);
    end
    if sum(contains(fieldnames(dBase),'cloud')) > 0
        cloud = [dBase.cloud]; % Percent
    else
        cloud = nan(1,nEns);
    end
    if ~SBA
        % Fill in with cloud index where empty
        whrNaN = isnan(cloud)';
        [wv,iwv] = find_nearest(750,wave);
        li750 = Li(:,iwv);
        es750 = Es(:,iwv);
        cloud(whrNaN & (li750./es750 < cloudIndexes(1))) = 0; % Clear
        cloud(whrNaN & ( li750./es750 >= cloudIndexes(1)) & li750./es750 < cloudIndexes(2) ) = 50; % Partly cloudy
        cloud(whrNaN & (li750./es750 > cloudIndexes(2))) = 100; % Fully overcast (no good for validation, but okay for AWR)
    end

    if sum(contains(fieldnames(dBase),'wind')) > 0
        wind = [dBase.wind];
    end
    if sum(contains(fieldnames(dBase),'SZA')) > 0
        sza = [dBase.SZA];
    else
        disp('Calculating SZA')        
        % sun_position is not vector-ready
        for i=1:length(dateTime)
            time.year = year(dateTime(i));
            time.month = month(dateTime(i));
            time.day = day(dateTime(i));
            time.hour = hour(dateTime(i));
            time.min = minute(dateTime(i));
            time.sec = second(dateTime(i));
            time.UTC = 0*second(dateTime(i));
            location.latitude = [dBase(i).latitude];
            location.longitude = [dBase(i).longitude];
        
            sun = sun_position(time, location);
            sza(i) = sun.zenith;
        end
    end
    if ~SBA
        if sum(contains(fieldnames(SASData.Ancillary),'REL_AZ')) > 0
            relAz = [dBase.REL_AZ];
        end
    end

    %% Set flags
    if validation
        flags.Cloud = (cloud > tCloud(1)) & (cloud < tCloud(2)); %Partly cloudy
    else
        flags.Cloud = false(1,length(cloud));
    end
    flags.Wind = wind > tWind;
    flags.SZA = sza < tSZA(1) | sza > tSZA(2);
    if ~SBA
        flags.RelAz = abs(relAz)<tRelAz(1) | abs(relAz)>tRelAz(2);
    end
    flags.QWIP = qwip > tQWIP;
    flags.QA = qa < tQA;
    waveRange = find(wave >= negRrs(1) & wave <= negRrs(2));
    flags.negRrs = any(Rrs(:,waveRange) < 0.0,2)';

    if ~SBA
        save(sprintf('dat/%s_flags.mat',cruise),"dateTime","wave","Rrs","Es","Li","Lw",...
            "chl","avw","qwip","qa","cloud","wind","sza","relAz","flags")
    else
        save(sprintf('dat/%s_flags.mat',cruise),"dateTime","wave","Rrs","Es","Lw",...
        "chl","avw","qwip","qa","cloud","wind","sza","flags")
    end


    %% Apply flags
    if ~SBA
        flag = [flags.Cloud] | [flags.Wind] | [flags.SZA] | [flags.RelAz] | [flags.QWIP] ...
            | [flags.negRrs];
    else
        % No RelAz flag
        flag = [flags.Cloud] | [flags.Wind] | [flags.SZA] | [flags.QWIP] ...
        | [flags.negRrs];
    end

    %% Write CSV file for awr2env.py
    % In effect, we will have 2 files. One will have 0 or 1, the other 0 or 2.
    % Round to the nearest second
    dateTime = dateshift(dateTime,'start','minute') + seconds(round(second(dateTime)));
    YEAR = year(dateTime)';
    MONTH = month(dateTime)';
    DAY = day(dateTime)';
    HOUR = hour(dateTime)';
    MINUTE = minute(dateTime)';
    SECOND = round(second(dateTime)');

    if validation
        csvOutFile = sprintf('dat/%s_flags.csv',cruise);
        FLAG = 2*int8(~flag'); % 0, 1, or 2 for reject, seabass-only, validation
    else
        csvOutFile = sprintf('dat/%s_all_flags.csv',cruise);
        FLAG = int8(~flag'); % 0, 1, or 2 for reject, seabass-only, validation
    end
    T = table(YEAR,MONTH,DAY,HOUR,MINUTE,SECOND,FLAG);
    writetable(T,csvOutFile)
else
    load(sprintf('dat/%s_flags.mat',cruise))
    nEns = length(Rrs); % check this
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

    if ~SBA
        ax3 = nexttile;
        plot(wave,Li,'k')
        hold on
        flagSpectra(ax3,wave,Li,flags,0)
        ylabel('L_i [\muW cm^{-2} nm^{-1} sr^{-1}]')
    end

    ax4 = nexttile;
    plot(wave,Lw,'k')
    hold on
    flagSpectra(ax4,wave,Lw,flags,0)
    ylabel('L_w [\muW cm^{-2} nm^{-1} sr^{-1}]')

    if all(isnan(cloud))
        th1 = text(0.70,0.9,sprintf('Cloud: Not reported'),'Units','normalized');
    else
        th1 = text(0.70,0.9,sprintf('Cloud: %d',sum(flags.Cloud)),'Units','normalized');
    end
    if all(isnan(wind))
        th2 = text(0.70,0.85,sprintf('Wind: Not reported'),'Units','normalized');
    else
        th2 = text(0.70,0.85,sprintf('Wind: %d',sum(flags.Wind)),'Units','normalized');
    end
    th3 = text(0.70,0.8,sprintf('SZA: %d',sum(flags.SZA)),'Units','normalized');
    if ~SBA
        th4 = text(0.70,0.75,sprintf('RelAz: %d',sum(flags.RelAz)),'Units','normalized');
    end
    th5 = text(0.70,0.7,sprintf('QWIP: %d',sum(flags.QWIP)),'Units','normalized');
    th6 = text(0.70,0.65,sprintf('QA: %d (not used)',sum(flags.QA)),'Units','normalized');%,'Color','r');
    th7 = text(0.70,0.60,sprintf('Neg. Rrs: %d',sum(flags.negRrs)),'Units','normalized');

    if ~SBA
        gud = ~flags.RelAz & ~flags.Wind & ~flags.SZA & ~flags.QWIP & ~flags.Cloud; %& ~flags.QA
        th8 = text(0.70,0.55,sprintf('Remaining: %d of %d',sum(gud),nEns),'Units','normalized');
        set([th1 th2 th3 th4 th5 th6 th7 th8],'FontName',fontName,'fontsize',12)
        set([ax1 ax2 ax3 ax4],'FontName',fontName,'FontSize',16, 'xgrid','on', 'ygrid','on')
    else
        gud = ~flags.Wind & ~flags.SZA & ~flags.QWIP & ~flags.Cloud; %& ~flags.QA
        th8 = text(0.70,0.55,sprintf('Remaining: %d of %d',sum(gud),nEns),'Units','normalized');
        set([th1 th2 th3 th5 th6 th7 th8],'FontName',fontName,'fontsize',12)
        set([ax1 ax2 ax4],'FontName',fontName,'FontSize',16, 'xgrid','on', 'ygrid','on')
    end
    set(fh3,'position',[1926         381        1196         979])

    if validation
        exportgraphics(fh3,sprintf('plt/%s_AllSpec_validation.png',cruise))
    else
        exportgraphics(fh3,sprintf('plt/%s_AllSpec.png',cruise))
    end

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
    hold on
    plot(dateTime([flags.QWIP]),qwip([flags.QWIP]),...
        'color','k','marker','x','markersize',18,...
        'linestyle','none');
    ylabel('QWIP')

    yyaxis(ax2,'right')
    ph4 = plot(dateTime,qa,'marker','.','markersize',18);
    ylabel('QA_score')

    ax3 = subplot(4,1,3);
    yyaxis(ax3,'left')
    ph5 = plot(dateTime,cloud,'marker','.','markersize',18);
    hold on
    plot(dateTime([flags.Cloud]),cloud([flags.Cloud]),...
        'color','k','marker','x','markersize',18,...
        'linestyle','none');
    ylabel('Cloud [%]')

    yyaxis(ax3,'right')
    ph6 = plot(dateTime,wind,'marker','.','markersize',18);
    hold on
    plot(dateTime([flags.Wind]),wind([flags.Wind]),...
        'color','k','marker','x','markersize',18,...
        'linestyle','none');
    ylabel('Wind [m/s]')

    ax4 = subplot(4,1,4);
    yyaxis(ax4,'left')
    ph7 = plot(dateTime,sza,'marker','.','markersize',18);
    hold on
    plot(dateTime([flags.SZA]),sza([flags.SZA]),...
        'color','k','marker','x','markersize',18,...
        'linestyle','none');
    ylabel('SZA')

    if ~SBA
        yyaxis(ax4,'right')
        ph8 = plot(dateTime,relAz,'marker','.','markersize',18);
        hold on
        plot(dateTime([flags.RelAz]),relAz([flags.RelAz]),...
            'color','k','marker','x','markersize',18,...
            'linestyle','none');
        ylabel('Relative Azimuth')
    end

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
        ah1 = axes;
        plot(wave,Rrs,'k')
        hold on
        flagSpectra(ah1,wave,Rrs,flags,0)
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

        if size(flag,1)~=1
            flags.Manual = flag';
        else
            flags.Manual = flag;
        end
        %% Apply flags
        if ~SBA
            flag = [flags.Cloud] | [flags.Wind] | [flags.SZA] | [flags.RelAz] | [flags.QWIP] |...
                [flags.Manual] | [flags.negRrs];
        else
            flag = [flags.Cloud] | [flags.Wind] | [flags.SZA] |  [flags.QWIP] |...
            [flags.Manual] | [flags.negRrs];
        end

        %% Write CSV file for awr2env.py
        % In effect, we will have 2 files. One will have 0 or 1, the other 0 or 2.
        % Round to the nearest second
        dateTime = dateshift(dateTime,'start','minute') + seconds(round(second(dateTime)));
        YEAR = year(dateTime)';
        MONTH = month(dateTime)';
        DAY = day(dateTime)';
        HOUR = hour(dateTime)';
        MINUTE = minute(dateTime)';
        SECOND = round(second(dateTime))';

        if validation
            csvOutFile = sprintf('dat/%s_flags.csv',cruise);
            FLAG = 2*int8(~flag'); % 0, 1, or 2 for reject, seabass-all, validation
        else
            csvOutFile = sprintf('dat/%s_all_flags.csv',cruise);
            FLAG = int8(~flag'); % 0, 1, or 2 for reject, seabass-all, validation
        end

        T = table(YEAR,MONTH,DAY,HOUR,MINUTE,SECOND,FLAG);
        writetable(T,csvOutFile)
    end
end

if ~SBA
    save(sprintf('dat/%s_flags.mat',cruise),"dateTime","wave","Rrs","Es","Li","Lw",...
        "chl","avw","qwip","qa","cloud","wind","sza","relAz","flags")
else
    save(sprintf('dat/%s_flags.mat',cruise),"dateTime","wave","Rrs","Es","Lw",...
        "chl","avw","qwip","qa","cloud","wind","sza","flags")
end


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
if sum(contains(fieldnames(flags),'RelAz')) > 0
    if sum(flags.RelAz)>0
        ph4 = plot(ax,wave,Var(flags.RelAz,:),'b','linewidth',2);
    else
        ph4 = plot(ax,wave,Var*nan,'b','linewidth',2);
    end
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
if sum(flags.negRrs)>0
    ph7 = plot(ax,wave,Var(flags.negRrs,:),'m','linewidth',2,'linestyle','--');
else
    ph7 = plot(ax,wave,Var*nan,'m','linewidth',2,'linestyle','--');
end


% These won't concatonate
% set([ph1 ph2 ph3 ph4 ph5 ph6],'linewidth',2)
if leg
    legend([ph1(1) ph2(1)  ph3(1) ph4(1)  ph5(1)  ph6(1) ph7(1)],...
        fieldnames(flags))
end

end

