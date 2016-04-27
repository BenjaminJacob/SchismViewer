function h=myPlotTr(st,x,z,xylims,clims)
% st - data mat to plot [np,nlevs]
% x - distance along the trasect
% z - levs [np,nlevs]

if min(size(st))==1
  h=plot(x,st);
else
  h=pcolor(x,z',st');
  set(h,'EdgeColor','none','FaceColor','interp');
end

if length(xylims)
  axis(xylims)
end
if length(clims)
  caxis(clims)
end

