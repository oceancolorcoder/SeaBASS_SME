% Input data from make_database_hypercp.m
% Add AVW and QWIP if not already included in Derived Products
% Display spectral (ir)radiance and reflectance with ancillary data
% dashboard
% Rank as   1) Validation, NOMAD, and SeaBASS
%           2) NOMAD and SeaBASS
%           3) SeaBASS only
%           4) Fail
%   These designations are subjective, but technically for 1 it would need
%   to be cloud free which would have to come from ancillary data.

% NOTE: Autonomous collections generally have many hundreds of spectra.
% There is no good, practical way to inspect every single one... This is
% going to be a problem.

% Screen on QWIP > 0.2
% Screen for cloud. Only validation if cloud=0%
% Screen on wind. <=7 m/s for validation (IOCCG 2019), <=10 m/s for NOMAD (NASA
% Protocols 2003), <=15 m/s SeaBASS (Zibordi et al. 2009)


wipe
fontName = machine_prefs;
%% Setup
cruise = 'EXPORTSNA_Boss';

%% Load and overview
load(sprintf('dat/%s.mat',cruise))
nEns = size([SASData.Ancillary.datenum],2);

Rrs = vertcat(SASData.Reflectance.Rrs_HYPER);
Es = vertcat(SASData.Irradiance.ES_HYPER);
Li = vertcat(SASData.Radiance.LI_HYPER);
Lw = vertcat(SASData.Radiance.LW_HYPER);

try
    % Potential of wavelength mismatch depending on glint processing
    wavelength = vertcat(SASData.Reflectance.Rrs_HYPER_wavelength);
catch
    disp('Wavelength cannot be concatonated')
    % In this case, interpolation of data to the shorter wavelength vector
    % will be required

end
wavelength = wavelength(1,:);

[abs_QWIP_score,QCI,AVW] = AVW_QWIP_2D_fun(Rrs,wavelength,'none','none');
minMaxPlot = [440 590];
fh12 = QWIP_figure_fun(Rrs,wavelength,AVW,cruise,minMaxPlot);

%% Screen on QWIP
Q0p2 = abs_QWIP_score > 0.2;

fh3 = figure;
th1 = tiledlayout(2,2);
ax1 = nexttile;
plot(wavelength,Rrs)
if any(Q0p2)
    hold
    plot(wavelength,Rrs(Q0p2,:),'Color','r','LineWidth',2)
end
ylabel('R_{rs} [sr^{-1}]')
set(gca,'FontName',fontName,'FontSize',16)
grid on
ax2 = nexttile;
plot(wavelength,Es)
if any(Q0p2)
    hold
    plot(wavelength,Es(Q0p2,:),'Color','r','LineWidth',2)
end
ylabel('E_s [\muW cm^{-2} nm^{-1}]')
set(gca,'FontName',fontName,'FontSize',16)
grid on
ax3 = nexttile;
plot(wavelength,Li)
if any(Q0p2)
    hold
    plot(wavelength,Li(Q0p2,:),'Color','r','LineWidth',2)
end
ylabel('L_i [\muW cm^{-2} nm^{-1} sr^{-1}]')
set(gca,'FontName',fontName,'FontSize',16)
grid on
ax4 = nexttile;
plot(wavelength,Lw)
if any(Q0p2)
    hold
    plot(wavelength,Lt(Q0p2,:),'Color','r','LineWidth',2)
end
ylabel('L_w [\muW cm^{-2} nm^{-1} sr^{-1}]')
set(gca,'FontName',fontName,'FontSize',16)
grid on

set(fh3,'position',[1926         381        1196         979])

exportgraphics(fh12(1),sprintf('plt/%s_QWIP.png',cruise))
exportgraphics(fh12(2),sprintf('plt/%s_QWIP_Hist.png',cruise))
exportgraphics(fh3,sprintf('plt/%s_AllSpec.png',cruise))



