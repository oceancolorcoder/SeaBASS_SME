function handles = plotFlags(AWR,ancillary,flags)

fontName = machine_prefs;

handles.fh3 = figure;
tile1 = tiledlayout(2,2);
handles.ax11 = nexttile;
plot(AWR.wave,AWR.Rrs,'k')
hold on
if exist('AWR.Rrs_sd','var')
    errorbar(AWR.wave_sd,AWR.Rrs_sub,AWR.Rrs_sd,'color','k')
end

flagSpectra(handles.ax11,AWR.wave,AWR.Rrs,flags,1)
ylabel('R_{rs} [sr^{-1}]')

handles.ax12 = nexttile;
plot(AWR.wave,AWR.Es,'k')
hold on
flagSpectra(handles.ax12,AWR.wave,AWR.Es,flags,0)
ylabel('E_s [\muW cm^{-2} nm^{-1}]')

if ~ancillary.SBA
    handles.ax13 = nexttile;
    plot(AWR.wave,AWR.Li,'k')
    hold on
    flagSpectra(handles.ax13,AWR.wave,AWR.Li,flags,0)
    ylabel('L_i [\muW cm^{-2} nm^{-1} sr^{-1}]')
end

handles.ax14 = nexttile;
plot(AWR.wave,AWR.Lw,'k')
hold on
flagSpectra(handles.ax14,AWR.wave,AWR.Lw,flags,0)
ylabel('L_w [\muW cm^{-2} nm^{-1} sr^{-1}]')

if all(isnan(ancillary.cloud))
    th1 = text(0.70,0.9,sprintf('Cloud: Not reported'),'Units','normalized');
else
    th1 = text(0.70,0.9,sprintf('Cloud: %d',sum(flags.Cloud)),'Units','normalized');
end
if all(isnan(ancillary.wind))
    th2 = text(0.70,0.85,sprintf('Wind: Not reported'),'Units','normalized');
else
    th2 = text(0.70,0.85,sprintf('Wind: %d',sum(flags.Wind)),'Units','normalized');
end
th3 = text(0.70,0.8,sprintf('SZA: %d',sum(flags.SZA)),'Units','normalized');
if ~ancillary.SBA
    th4 = text(0.70,0.75,sprintf('RelAz: %d',sum(flags.RelAz)),'Units','normalized');
end
th5 = text(0.70,0.7,sprintf('QWIP: %d',sum(flags.QWIP)),'Units','normalized');
th6 = text(0.70,0.65,sprintf('QA: %d (not used)',sum(flags.QA)),'Units','normalized');%,'Color','r');
th7 = text(0.70,0.60,sprintf('Neg. Rrs: %d',sum(flags.negRrs)),'Units','normalized');

if ~ancillary.SBA
    gud = ~flags.RelAz & ~flags.Wind & ~flags.SZA & ~flags.QWIP & ~flags.Cloud; %& ~flags.QA
    handles.th8 = text(0.70,0.55,sprintf('Remaining: %d of %d',sum(gud),AWR.nSpectra),'Units','normalized');
    set([th1 th2 th3 th4 th5 th6 th7 handles.th8],'FontName',fontName,'fontsize',12)
    set([ax11 ax12 ax13 ax14],'FontName',fontName,'FontSize',16, 'xgrid','on', 'ygrid','on')
else
    gud = ~flags.Wind & ~flags.SZA & ~flags.QWIP & ~flags.Cloud; %& ~flags.QA
    handles.th8 = text(0.70,0.55,sprintf('Remaining: %d of %d',sum(gud),AWR.nSpectra),'Units','normalized');
    set([th1 th2 th3 th5 th6 th7 handles.th8],'FontName',fontName,'fontsize',12)
    set([handles.ax11 handles.ax12 handles.ax14],'FontName',fontName,'FontSize',16, 'xgrid','on', 'ygrid','on')
end
set(handles.fh3,'position',[1926         381        1196         979])

if ancillary.validation
    exportgraphics(handles.fh3,sprintf('plt/%s_spec.png',ancillary.cruise))
else
    exportgraphics(handles.fh3,sprintf('plt/%s_all_spec.png',ancillary.cruise))
end

%% Figure timeline Chl, AVW, QWIP, Wei score, ..
fh4 = figure;
set(fh4,'Position',[200 200 850 950])
ax21 = subplot(4,1,1);
title(strrep(ancillary.cruise,'_','-'))
yyaxis(ax21,'left')
ph1 = plot(AWR.dateTime,ancillary.chl,'marker','.','markersize',18);
ylabel('chlor_a')

yyaxis(ax21,'right')
ph2 = plot(AWR.dateTime,ancillary.avw,'marker','.','markersize',18);
ylabel('AVW [nm]')

ax22 = subplot(4,1,2);
yyaxis(ax22,'left')
ph3 = plot(AWR.dateTime,ancillary.qwip,'marker','.','markersize',18);
hold on
plot(AWR.dateTime([flags.QWIP]),ancillary.qwip([flags.QWIP]),...
    'color','k','marker','x','markersize',18,...
    'linestyle','none');
ylabel('QWIP')

yyaxis(ax22,'right')
ph4 = plot(AWR.dateTime,ancillary.qa,'marker','.','markersize',18);
ylabel('QA_score')

ax23 = subplot(4,1,3);
yyaxis(ax23,'left')
ph5 = plot(AWR.dateTime,ancillary.cloud,'marker','.','markersize',18);
hold on
plot(AWR.dateTime([flags.Cloud]),ancillary.cloud([flags.Cloud]),...
    'color','k','marker','x','markersize',18,...
    'linestyle','none');
ylabel('Cloud [%]')

yyaxis(ax23,'right')
ph6 = plot(AWR.dateTime,ancillary.wind,'marker','.','markersize',18);
hold on
plot(AWR.dateTime([flags.Wind]),ancillary.wind([flags.Wind]),...
    'color','k','marker','x','markersize',18,...
    'linestyle','none');
ylabel('Wind [m/s]')

ax24 = subplot(4,1,4);
yyaxis(ax24,'left')
ph7 = plot(AWR.dateTime,ancillary.sza,'marker','.','markersize',18);
hold on
plot(AWR.dateTime([flags.SZA]),ancillary.sza([flags.SZA]),...
    'color','k','marker','x','markersize',18,...
    'linestyle','none');
ylabel('SZA')

if ~ancillary.SBA
    yyaxis(ax24,'right')
    ph8 = plot(AWR.dateTime,ancillary.relAz,'marker','.','markersize',18);
    hold on
    plot(AWR.dateTime([flags.RelAz]),ancillary.relAz([flags.RelAz]),...
        'color','k','marker','x','markersize',18,...
        'linestyle','none');
    ylabel('Relative Azimuth')
end

set([ax21 ax22 ax23 ax24],'fontname',fontName, 'fontsize', 14, 'xgrid', 'on')
if ancillary.validation
    exportgraphics(fh4,sprintf('plt/%s_timeline.png',ancillary.cruise))
else
    exportgraphics(fh4,sprintf('plt/%s_all_timeline.png',ancillary.cruise))
end