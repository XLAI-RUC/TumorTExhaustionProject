function [x1_n,x2_n,y1_n,y2_n] = normalize_f(ax1,axPos,x1,x2,y1,y2)
x1_n = axPos(1) + (x1 - ax1.XLim(1)) / diff(ax1.XLim) * axPos(3);
y1_n = axPos(2) + (y1- ax1.YLim(1)) / diff(ax1.YLim) * axPos(4);
x2_n = axPos(1) + (x2 - ax1.XLim(1)) / diff(ax1.XLim) * axPos(3);
y2_n = axPos(2) + (y2- ax1.YLim(1)) / diff(ax1.YLim) * axPos(4);
end