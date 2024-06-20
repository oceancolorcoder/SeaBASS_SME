function manualFlag(ancillary,handles,AWR,flags)

if ancillary.validation
    fprintf('Screen for validation\n')
else
    fprintf('Screen for NOMAD\n')
end

flag = zeros(1,AWR.nSpectra);

disp('Manual Screening')
disp('Zoom to spectrum of interest and hit Continue')
disp('Left mouse to flag spectrum. Accept to save and move on. Reject to ignore last.')
disp('Middle mouse to exit and save.')

button1 = uicontrol(handles.fh5,'Style', 'pushbutton', 'String', 'Continue after zooming',...
    'Position', [5 150 130 25], 'Callback', @buttonCallBack);
button2 = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
    'Position', [5 100 50 25], 'Callback', @buttonCallBack);
button3 = uicontrol('Style', 'pushbutton', 'String', 'Reject',...
    'Position', [5 50 50 25], 'Callback', @buttonCallBack);
button4 = uicontrol('Style', 'pushbutton', 'String', 'Exit',...
    'Position', [5 5 50 25], 'Callback', @buttonCallBack);

% for tries=1:10
while 1    
    disp('Zoom to spectrum of interest and hit Continue')
    
    uiwait(handles.fh5)
    if get(button4,'userdata') == 1
        disp('Exitting and saving')
        break
    end

    [x,y] = ginput(1);
    plot(x,y,'k*')

    set(button1,'userdata',0)
    set(button2,'userdata',0)
    set(button3,'userdata',0)
    set(button4,'userdata',0)
    
    uiwait

    if get(button1,'userdata') == 1 || get(button3,'userdata') == 1
        set(button1,'userdata',0)
        set(button2,'userdata',0)
        set(button3,'userdata',0)
        set(button4,'userdata',0)
        continue
    elseif get(button2,'userdata') == 1
        [~,windex] = find_nearest(x,AWR.wave);
        RrsX = AWR.Rrs(:,windex);
        [~,Rindex] = find_nearest(y,RrsX);
        plot(AWR.wave,AWR.Rrs(Rindex,:),'k','LineWidth',3)
        flag(Rindex) = 1;
        fprintf('Index of selected spectrum: %d\n',Rindex)
        % disp(flag(Rindex))
        set(button1,'userdata',0)
        set(button2,'userdata',0)
        set(button3,'userdata',0)
        set(button4,'userdata',0)
    elseif get(button4,'userdata') == 1
        disp('Exitting and saving')
        break
    end
end

if size(flag,1)~=1
    flags.Manual = logical(flag');
else
    flags.Manual = logical(flag);
end
%% Apply flags
if ~ancillary.SBA
    flag = [flags.Cloud] | [flags.Wind] | [flags.SZA] | [flags.RelAz] | [flags.QWIP] |...
        [flags.Manual] | [flags.negRrs];
    set(handles.th8,'String',sprintf('Remaining: %d of %d',sum(~flag),AWR.nSpectra));
else
    flag = [flags.Cloud] | [flags.Wind] | [flags.SZA] |  [flags.QWIP] |...
        [flags.Manual] | [flags.negRrs];
    set(handles.th8,'String',sprintf('Remaining: %d of %d',sum(~flag),AWR.nSpectra));
end

flagSpectra(handles.ax11,AWR.wave,AWR.Rrs,flags,1)
flagSpectra(handles.ax12,AWR.wave,AWR.Es,flags,0)
if ~ancillary.SBA
    flagSpectra(handles.ax13,AWR.wave,AWR.Li,flags,0)
end

if sum(contains(fieldnames(AWR),'Lw')) > 0
    flagSpectra(handles.ax14,AWR.wave,AWR.Lw,flags,0)
elseif sum(contains(fieldnames(AWR),'Lt')) > 0
    flagSpectra(handles.ax14,AWR.wave,AWR.Lt,flags,0)
end

if ancillary.validation
    exportgraphics(handles.fh3,sprintf('plt/%s_spec.png',ancillary.cruise))
else
    exportgraphics(handles.fh3,sprintf('plt/%s_all_spec.png',ancillary.cruise))
end

%% Write CSV file for awr2env.py
% In effect, we will have 2 files. One will have 0 or 1, the other 0 or 2.
% Round to the nearest second
dateTime = dateshift(AWR.dateTime,'start','minute') + seconds(round(second(AWR.dateTime)));
[Yr,Mon,Day,Hr,Min,Sec] = datevec(dateTime');

if ancillary.validation
    csvOutFile = sprintf('dat/%s_flags.csv',ancillary.cruise);
    FLAG = 2*int8(~flag'); % 0, 1, or 2 for reject, seabass-only, validation
else
    csvOutFile = sprintf('dat/%s_all_flags.csv',ancillary.cruise);
    FLAG = int8(~flag'); % 0, 1, or 2 for reject, seabass-only, validation
end
T = table(Yr,Mon,Day,Hr,Min,Sec,FLAG);
writetable(T,csvOutFile)

end

function buttonCallBack(hObject,eventdata)
set(hObject,'UserData',1)
uiresume
end
