function testButton

if ancillary.validation
    fprintf('Screen for validation\n')
else
    fprintf('Screen for NOMAD\n')
end

% input('Continue now? (enter)');
% close all
fh5 = figure;
set(fh5,'position',[1926         381        1196         979])
ah1 = axes;
plot(AWR.wave,AWR.Rrs,'k')
hold on
flagSpectra(ah1,AWR.wave,AWR.Rrs,flags,0)
ylabel('R_{rs} [sr^{-1}]')

disp('Manual Screening')
disp('Zoom to spectrum of interest and hit Continue')
disp('Left mouse to flag spectrum. Continue to save and move on. Try again to ignore last.')
disp('Middle mouse to exit and save.')
flag = zeros(1,AWR.nSpectra);
for tries=1:10
    % h1 = uicontrol(fh5,'Style', 'pushbutton', 'String', 'Continue',...
    %     'Position', [5 100 50 25], 'Callback', 'butt=1; uiresume');
    h1 = uicontrol(fh5,'Style', 'togglebutton', 'String', 'Continue',...
        'Position', [5 100 50 25], 'Callback', 'butt=1; uiresume');

    % 
    uicontrol(h1)
    uiwait
    disp(butt)