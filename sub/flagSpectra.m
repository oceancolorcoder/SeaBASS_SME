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
% if sum(flags.QA)>0
%     ph6 = plot(ax,wave,Var(flags.QA,:),'b','linewidth',2,'linestyle','--');
% else
%     ph6 = plot(ax,wave,Var*nan,'b','linewidth',2,'linestyle','--');
% end
if sum(flags.negRrs)>0
    ph7 = plot(ax,wave,Var(flags.negRrs,:),'color','b','linewidth',2,'linestyle','--');
else
    ph7 = plot(ax,wave,Var*nan,'color','b','linewidth',2,'linestyle','--');
end
if sum(flags.Manual)>0
    ph8 = plot(ax,wave,Var(flags.Manual,:),'color','m','linewidth',2,'linestyle','--');
else
    ph8 = plot(ax,wave,Var*nan,'color','m','linewidth',2,'linestyle','--');
end


% These won't concatonate
% set([ph1 ph2 ph3 ph4 ph5 ph6],'linewidth',2)
if leg
    legNames = fieldnames(flags);
    if sum(contains(fieldnames(flags),'RelAz')) > 0        
        legend([ph1(1) ph2(1)  ph3(1) ph4(1)  ph5(1)  ph7(1) ph8(1)],...
            legNames([1:5 7 8]))
    else
        legend([ph1(1) ph2(1)  ph3(1) ph5(1)  ph7(1) ph8(1)],...
            legNames([1:4 6 7]))
    end
end

end