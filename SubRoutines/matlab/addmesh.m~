function ph_mesh=addmesh(gr,varargin)

axe=axis;

% dx=0;
% xmin=(axe(1))-dx;
% xmax=(axe(2))+dx;
% ymin=(axe(3))-dx;
% ymax=(axe(4))+dx;

xlims=xlim;
xmin=xlims(1);
xmax=xlims(2);

ylims=ylim;
ymin=ylims(1);
ymax=ylims(2);





% in nodes
innodes=find(gr.hgrid.nodes(:,2)  > xmin & xmax > gr.hgrid.nodes(:,2) &...
gr.hgrid.nodes(:,3)  > ymin & ymax > gr.hgrid.nodes(:,3)); 


area_cls=find(sum(ismember(gr.hgrid.elem(:,3:6),innodes),2));

triangles=gr.hgrid.elem(area_cls(:),2)==3;
quads=gr.hgrid.elem(area_cls(:),2)==4;


hold on
p = patch('Faces',gr.hgrid.elem(area_cls(quads),3:6),'Vertices',[gr.hgrid.x,gr.hgrid.y],...
    'FaceColor','none','EdgeColor','m'); %Plot Uqds
ph_mesh=triplot(gr.hgrid.elem(area_cls(triangles),3:5),gr.hgrid.x,gr.hgrid.y,'k'); %Plot Triangles
hold off


%add numbers
if nargin==2
    
    hold on
    
    %triangles
    nodes=unique(gr.hgrid.elem(area_cls,3:5));
    text(gr.hgrid.x(nodes),gr.hgrid.y(nodes),num2str(nodes))
    xx=mean(gr.hgrid.x(gr.hgrid.elem(area_cls(triangles),3:5)),2);
    yy=mean(gr.hgrid.y(gr.hgrid.elem(area_cls(triangles),3:5)),2);
    plot(xx,yy,'r+')
    text(xx,yy,num2str(area_cls(triangles)),'color','r')
    
    
    %quads
    nodes=unique(gr.hgrid.elem(area_cls(quads),3:6));
    text(gr.hgrid.x(nodes),gr.hgrid.y(nodes),num2str(nodes))
    xx=mean(gr.hgrid.x(gr.hgrid.elem(area_cls(quads),3:6)),2);
    yy=mean(gr.hgrid.y(gr.hgrid.elem(area_cls(quads),3:6)),2);
    plot(xx,yy,'r+')
    text(xx,yy,num2str(area_cls(quads)),'color','r')
end




end