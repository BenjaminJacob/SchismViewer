function []=gr_plot2(gr,c,maxmin)
%function [h]=gr_plot2(gr,c,[maxmin])
% displays individula grid of gr=readGrid structure
% gr is one of the gr.hgrid, gr.ecenters, ....
% c - index color matrix for each node
% maxmin[max min] -maximum and minimum of the color stretch
% Example: gr_plot2(gr.ecenters,gr.ecenters.depth)
%
% Sergey Frolov,   March 2004

r=max(c)-min(c);
ncol=200;

if nargin ==3
   c( find(c < maxmin(1)) )= maxmin(1);
   c( find(c > maxmin(2)) )= maxmin(2);   
end
% c=histeq(c/r,ncol);

close
h=figure;

patch('Faces',gr.elem(find(gr.elem(:,2)==3),3:end-1),'Vertices',gr.nodes(:,2:end-1),'FaceVertexCData',c,'FaceColor','interp','EdgeColor','none')
patch('Faces',gr.elem(find(gr.elem(:,2)==4),3:end),'Vertices',gr.nodes(:,2:end-1),'FaceVertexCData',c,'FaceColor','interp','EdgeColor','none')
set(gca,'FontSize',14);

if nargin ==3
	caxis(maxmin)
end

colormap(jet(ncol))
h=colorbar;
set(h,'FontSize',14);
axis equal
