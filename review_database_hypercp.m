wipe

%% Setup
cruise = 'viirs_2019_foster'; % Has AVW but not QWIP yet
[fontName,projPath,imgPath] = machine_prefs();

inDbase = sprintf('dat/%s.mat',cruise);

%%
load(inDbase)
dateTime = datetime([SASData.Ancillary.datenum],'ConvertFrom','datenum','TimeZone','UTC');

%% Figure 1; timeline Chl, AVW, QWIP, Wei score, ..

fh1 = figure;
set(fh1,'Position',[200 200 850 950])
ax1 = subplot(4,1,1);
ph1 = plot(dateTime,[SASData.Derived_Products.chlor_a],'marker','.','markersize',18);
ylabel('chlor_a')
set(ax1,'fontname',fontName, 'fontsize', 14)

if sum(contains(fieldnames(SASData.Derived_Products),'avw')) > 0
    ax2 = subplot(4,1,2);
    ph2 = plot(dateTime,[SASData.Derived_Products.avw],'marker','.','markersize',18);
    ylabel('AVW [nm]')
    set(ax2,'fontname',fontName, 'fontsize', 14)
end

if sum(contains(fieldnames(SASData.Derived_Products),'qwip')) > 0
    ax3 = subplot(4,1,3);
    ph3 = plot(dateTime,[SASData.Derived_Products.qwip],'marker','.','markersize',18);
    ylabel('QWIP')
    set(ax3,'fontname',fontName, 'fontsize', 14)
else
    % Calculate QWIP
end

if sum(contains(fieldnames(SASData.Ancillary),'CLOUD')) > 0
    ax4 = subplot(4,1,4);
    ph4 = plot(dateTime,[SASData.Ancillary.CLOUD],'marker','.','markersize',18);
    ylabel('Cloud [%]')
    set(ax4,'fontname',fontName, 'fontsize', 14)
else
    % test cloud indexes
end
