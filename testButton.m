
wipe
f = figure('Position',[500 500 400 300]);
data.value = 0;
fprintf('value = %d\n',data.value)
guidata(f,data)
testButt(f)
% guidata(f,data)
fprintf('value = %d\n',data.value)

function testButt(src,event)

uicontrol('String','Continue','Callback','uiresume(f)');
uiwait(src)
data = guidata(src);
fprintf('value = %d\n',data.value)
data.value = 1;
fprintf('value = %d\n',data.value)
guidata(src,data)
% data
end
