function [h]=gr_plot(gr,maxmin)
%function [h]=gr_plot(gr,[maxmin])
% displays a grid gr.hgrid
% requires gr.hgrid.nodes, gr.hgrid.depth , gr.hgrid.elem
% maxmin[max min] -maximum and minimum of the color stretch
%
% Sergey Frolov,   March 2004

I=gr.hgrid.depth(:);
r=max(I)-min(I);
ncol=200;

if nargin ==2
   I( find(I < maxmin(1)) )= maxmin(1);
   I( find(I > maxmin(2)) )= maxmin(2);   
end
% I=histeq(I/r,ncol);

close
h=figure;

patch('Faces',gr.hgrid.elem(find(gr.hgrid.elem(:,2)==3),3:end-1),'Vertices',gr.hgrid.nodes(:,2:end-1),'FaceVertexCData',I,'FaceColor','interp','EdgeColor','b')
patch('Faces',gr.hgrid.elem(find(gr.hgrid.elem(:,2)==4),3:end),'Vertices',gr.hgrid.nodes(:,2:end-1),'FaceVertexCData',I,'FaceColor','interp','EdgeColor','b')

colormap(jet(ncol))
colorbar
axis equal