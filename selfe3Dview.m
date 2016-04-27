

%% SELFE 3 vie Plot


buro='/media/benjamin/TOSHIBA EXT/Arbeit/';
hause='/media/benjamin/Pladde/';
ort=buro;

fn.runDir=[ort 'Promotion/Selfe_Setup_reference/'];   %grid verzeichniss
fn.outDir=[ort 'Promotion/Selfe_Setup_reference/combined/'];   %Lage der Outpudaten
fn.TideAna='/home/benjamin/Desktop/SELFE_Modellierung/Bathymetry_experimente/scenarioA/TideAna'; %Lage derharmonischen Analyse
%%

grName='hgrid.gr3';
gr.hgrid=gr_readHGrid( fullfile(fn.runDir,grName) );


gr.hgrid0=gr.hgrid;
gr.hgrid=ll.hgrid;

gr_plot(gr.hgrid,gr.hgrid.depth)


axe=axis;

xlims=xlim;
xmin=xlims(1);
xmax=xlims(2);
ylims=ylim;
ymin=ylims(1);
ymax=ylims(2);



 area_cls=find(sum(gr.hgrid.x(gr.hgrid.elem(:,3:5)) > xmin & xmax > gr.hgrid.x(gr.hgrid.elem(:,3:5)),2)>=1 ...
      & sum(gr.hgrid.y(gr.hgrid.elem(:,3:5)) > ymin & ymax > gr.hgrid.y(gr.hgrid.elem(:,3:5)),2)>=1 );

  
 trinagles=gr.hgrid.elem(area_cls(:),2)==3;
 quads=gr.hgrid.elem(area_cls(:),2)==4;
 
 figure
 hold on
 ph_mesh=triplot(gr.hgrid.elem(area_cls(trinagles),3:5),gr.hgrid.x,gr.hgrid.y,'k'); %Plot Triangles
 
 
 
stack=13;

%elevation
binaryfile='elev.61';
file=sprintf('%i_%s',stack,binaryfile);
file=fullfile(fn.outDir,file);
helev=sz_readHeader(file);

% %
% binaryfile='zcor.63';
% file=sprintf('%i_%s',stack,binaryfile);
% file=fullfile(fn.outDir,file);
% hzcor=sz_readHeader(file);

%hvel 
binaryfile='hvel.64';
file=sprintf('%i_%s',stack,binaryfile);
file=fullfile(fn.outDir,file);
hhvel=sz_readHeader(file);


%hvel 
binaryfile='vert.63';
file=sprintf('%i_%s',stack,binaryfile);
file=fullfile(fn.outDir,file);
hvert=sz_readHeader(file);


timesteps=1:helev.nSteps;

[elev ts]=sz_readTimeStep(helev,timesteps);


%[zcor ts]=sz_readTimeStep(hzcor,timesteps);

[hvel ts]=sz_readTimeStep(hhvel,timesteps);


[vert ts]=sz_readTimeStep(hvert,timesteps);

ti=24


   %axis: 7.5060    8.0532   53.6489   54.0122   -0.3785    0.0541
 figure
 trisurf(gr.hgrid.elem(area_cls(trinagles),3:5),gr.hgrid.x,gr.hgrid.y,-gr.hgrid.depth/100);
 shading interp
 caxis([-25 3]/100)
 %colormap copper
%shading interp
colormap haxby
colorbar 

%hlight=light('Position',[1 1 0],'Style','infinite');
hlight=light('Position',[7.3021   54.0068    2.6008],'Style','infinite');

%cpos=campos

minx=min(min(gr.hgrid.x(gr.hgrid.elem(area_cls(trinagles),3:5))));
maxx=max(max(gr.hgrid.x(gr.hgrid.elem(area_cls(trinagles),3:5))));
miny=min(min(gr.hgrid.y(gr.hgrid.elem(area_cls(trinagles),3:5))));
maxy=max(max(gr.hgrid.y(gr.hgrid.elem(area_cls(trinagles),3:5))));
maxz=max(max(gr.hgrid.depth(gr.hgrid.elem(area_cls(trinagles),3:5))));
minz=min(min(gr.hgrid.depth(gr.hgrid.elem(area_cls(trinagles),3:5))));

dn=30;
dz=5;


dx=0.01;
xinterp=minx:dx:maxx;
dy=0.01;
yinterp=miny:dy:maxy;
zinterp=-(minz:3:maxz);

[XX YY ZZ]=meshgrid(xinterp,yinterp,zinterp);

%views=view


axis vis3d
grid off
freezeColors
hold on
for ti=timesteps
%set(hSSE,'EdgeColor','none','FAceColor',[0.2 0.2 0.8],'FaceAlpha',0.5)
title(['t = ' num2str(ti)])
zlevels=sz_computeZlevels(gr.hgrid.depth,elev(:,1,ti),helev.vgrid);

hSSE=trisurf(gr.hgrid.elem(area_cls(trinagles),3:5),gr.hgrid.x,gr.hgrid.y,elev(:,1,ti)/100,...
    'EdgeColor','none','FAceColor',[0.2 0.2 0.8],'FaceAlpha',0.6);

in_nodes=unique(gr.hgrid.elem(area_cls(trinagles),3:5));

 u=map_sz2hts(hhvel,hvel(:,1,ti),1);%
 v=map_sz2hts(hhvel,hvel(:,2,ti),1);
 w=map_sz2hts(hvert,vert(:,1,ti),1);

 
 
 
xx=repmat(gr.hgrid.x(in_nodes),1,helev.vgrid.nLevels);
yy=repmat(gr.hgrid.y(in_nodes),1,helev.vgrid.nLevels);
zz=zlevels(in_nodes,:);
 


uu=u(in_nodes,:)/110; 
vv=v(in_nodes,:)/110;
ww=w(in_nodes,:);%

Fu = TriScatteredInterp(xx(:),yy(:),zz(:),uu(:));
Fv = TriScatteredInterp(xx(:),yy(:),zz(:),vv(:));
Fw = TriScatteredInterp(xx(:),yy(:),zz(:),ww(:));


UU=Fu(XX,YY,ZZ);
VV=Fv(XX,YY,ZZ);
WW=Fw(XX,YY,ZZ);

fac=1;
zfac=0.01;
hcones=coneplot(XX*fac,YY*fac,ZZ*zfac,UU,VV,WW,XX*fac,YY*fac,ZZ*zfac,0.5);
set(hcones,'FaceColor','r','EdgeColor','r')

%figure
%quiver3(XX,YY,ZZ,UU,VV,WW);

%ww=zeros(size(uu));


%slice(XX,YY,ZZ*zfac,UU,[7.5:.25:8] ,[],[])
%colormap jet

%figure
%coneplot(xx(1:dn:end,1:dz:end),yy(1:dn:end,1:dz:end),zz(1:dn:end,1:dz:end),uu(1:dn:end,1:dz:end),vv(1:dn:end,1:dz:end),ww(1:dn:end,1:dz:end));
%hq=quiver3(xx(1:dn:end,1:dz:end),yy(1:dn:end,1:dz:end),zz(1:dn:end,1:dz:end),...
 %   uu(1:dn:end,1:dz:end),vv(1:dn:end,1:dz:end),ww(1:dn:end,1:dz:end),0,'color','r');

%flight
x0=7.5;
y0=54.2;

for phi=0:pi/5:2*pi
    dx=cos(phi)*2;
    dy=sin(phi)*2;
    campos([x0+dx y0+dy 2.5]);
    pause(1);
end
 
pause(0.4)
delete(hSSE)
%delete(hq)
delete(hcones)
end

zvec=-23:3:5;
[xx2 yy2 zz2]=meshgrid(sort(gr.hgrid.x(in_nodes(1:dn:end))),sort(gr.hgrid.y(in_nodes(1:dn:end))),zvec);

figure
coneplot(xx2,yy2,zz2,xx2,yy2,zz2,Fu(xx2,yy2,zz2),Fv(xx2,yy2,zz2),Fw(xx2,yy2,zz2))

figure
streamline(xx2,yy2,zz2,xx2,yy2,zz2,Fu(xx2,yy2,zz2),Fv(xx2,yy2,zz2),Fw(xx2,yy2,zz2))
 


figure
coneplot(XX,YY,ZZ,UU,VV,WW,XX,YY,ZZ)

figure
streamribbon(XX,YY,ZZ,UU,VV,WW,XX,YY,ZZ)

zslice=0;
xslice=mean(XX(1,:,1));
yslice=mean(YY(:,1,1));
zslice=mean(ZZ(1,1,:));

figure
slice(XX,YY,ZZ,UU,xslice,yslice,zslice)