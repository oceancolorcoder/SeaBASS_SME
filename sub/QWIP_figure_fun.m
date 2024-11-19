function fhs = QWIP_figure_fun(Rrs,wavelengths,AVW,Title,minMaxPlot)
% Adapted from R. Vandermuelen by D.Aurin 2022-09-28
% NDI = Normalized Difference Index = QCI (Quality Control Index?)


fontName = machine_prefs;

% Ensure vector in rows
AVW = AVW(:);

% Clean up title
Title = strrep(Title,'_','');

% minMaxPlot = [440 590]; % For plotting purposes

minMaxWave = [400 700]; % For truncating Rrs
[~,min_index] = min(abs(wavelengths-minMaxWave(1)));
[~,max_index] = min(abs(wavelengths-minMaxWave(2)));

RRS = Rrs(:,min_index:max_index);
wave_vis = wavelengths(:,min_index:max_index);

[~,index_490] = min(abs(wave_vis-490));
[~,index_665] = min(abs(wave_vis-665));
avw_bins = minMaxPlot(1):minMaxPlot(2);
QCI = (RRS(:,index_665) - RRS(:,index_490))./(RRS(:,index_665) + RRS(:,index_490));
p = [-8.399885e-09,1.715532e-05,-1.301670e-02,4.357838,-5.449532e02];
fit1 = polyval(p,avw_bins);

%% QCI/NDI v AVW & QWIP
fit1a = fit1 + 0.1;
fit1b = fit1 - 0.1;
fit2a = fit1 + 0.2;
fit2b = fit1 - 0.2;
fit3a = fit1 + 0.3;
fit3b = fit1 - 0.3;
fit4a = fit1 + 0.4;
fit4b = fit1 - 0.4;

fh1 = figure;

% Data
ph1 = plot(AVW,QCI,'LineStyle','none','Marker','o','MarkerSize',7,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5 0.5 1.0]);
hold on;

% Index lines
plot(avw_bins,fit1,'-k','LineWidth',2);
lh1 = plot(avw_bins,fit1a,'--g','LineWidth',2);
plot(avw_bins,fit1b,'--g','LineWidth',2);
lh2 = plot(avw_bins,fit2a,'--','LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
plot(avw_bins,fit2b,'--','LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
lh3 = plot(avw_bins,fit3a,'--','LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
plot(avw_bins,fit3b,'--','LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
lh4 = plot(avw_bins,fit4a,'-r','LineWidth',2);
plot(avw_bins,fit4b,'-r','LineWidth',2);
xlabel('AVW (nm)','FontSize',16);
ylabel(['NDI (' num2str(wave_vis(index_490)) ',' num2str(wave_vis(index_665)) ')'],'FontSize',16);
ylim([-2.0 1.5]);
xlim([minMaxPlot(1) minMaxPlot(2)]);

L{1} = 'QWIP ± 0.1';
L{2} = 'QWIP ± 0.2';
L{3} = 'QWIP ± 0.3';
L{4} = 'QWIP ± 0.4';
legh=legend([lh1 lh2 lh3 lh4], L,'Location','northwest','FontSize',12);
title(Title, 'FontSize',18)
text(550, 1.35,sprintf('n = %d',length(QCI)), 'FontSize',14)
set(gca,'FontName',fontName)


%% Histograms
% Make incremental bins of QWIP scores 
%   (0 - 0.1, 0.1 - 0.2, 0.2 - 0.3, 0.3 - 0.4, and > 0.4).
%   Then make some histograms to see the frequency of the QWIP scores
%   as a function of AVW. 

n = length(avw_bins);
hist_bin1 = zeros(n,1);
hist_bin2 = zeros(n,1);
hist_bin3 = zeros(n,1);
hist_bin4 = zeros(n,1);
hist_bin5 = zeros(n,1);

start = minMaxPlot(1) -1;
% for i = minMaxPlot(1):1:minMaxPlot(2)
% start = 399;
for i = avw_bins
    j = i+0.9999;    
    flag = AVW>i & AVW<j;

    flag2 = QCI>=fit1(i-start)-0.1 & QCI<=fit1(i-start)+0.1;
    flag3 = flag.*flag2;
    hist_bin1(i-start,1)=sum(flag3,'all');

    flag2 = (QCI>=fit1(i-start)-0.2 & QCI<fit1(i-start)-0.1) + (QCI<=fit1(i-start)+0.2 & QCI>fit1(i-start)+0.1);
    flag3 = flag.*flag2;
    hist_bin2(i-start,1)=sum(flag3,'all');

    flag2 = (QCI>=fit1(i-start)-0.3 & QCI<fit1(i-start)-0.2) + (QCI<=fit1(i-start)+0.3 & QCI>fit1(i-start)+0.2);
    flag3 = flag.*flag2;
    hist_bin3(i-start,1)=sum(flag3,'all');

    flag2 = (QCI>=fit1(i-start)-0.4 & QCI<fit1(i-start)-0.3) + (QCI<=fit1(i-start)+0.4 & QCI>fit1(i-start)+0.3);
    flag3 = flag.*flag2;
    hist_bin4(i-start,1)=sum(flag3,'all');

    flag2 = (QCI<fit1(i-start)-0.4) + (QCI>fit1(i-start)+0.4);
    flag3 = flag.*flag2;
    hist_bin5(i-start,1)=sum(flag3,'all');

end

fh2 = figure;
scaleFac = max(max([hist_bin1 hist_bin2 hist_bin3 hist_bin4 hist_bin5]));
bar((avw_bins),hist_bin1/scaleFac,'g');hold on;
bar((avw_bins),hist_bin2/scaleFac,'FaceColor',[0.9290 0.6940 0.1250])
bar((avw_bins),hist_bin3/scaleFac,'FaceColor',[0.8500 0.3250 0.0980])
bar((avw_bins),hist_bin4/scaleFac,'r')
bar((avw_bins),hist_bin5/scaleFac,'FaceColor',[0.6350 0.0780 0.1840])
xlabel('AVW (nm)','FontSize',16);
ylabel(sprintf('Frequency x %d',scaleFac),'FontSize',16);
legend(['QWIP = |0.0 - 0.1|, n = ' num2str(sum(hist_bin1))],['QWIP = |0.1 - 0.2|, n = ' num2str(sum(hist_bin2))],['QWIP = |0.2 - 0.3|, n = ' num2str(sum(hist_bin3))],['QWIP = |0.3 - 0.4|, n = ' num2str(sum(hist_bin4))],['QWIP > |0.4|, n = ' num2str(sum(hist_bin5))],'Location','northeast');
grid on
xlim([minMaxPlot(1) minMaxPlot(2)]);
hold off
title(Title, 'FontSize',18)
set(gca,'fontname',fontName)

fhs = [fh1 fh2];