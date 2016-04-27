
function transect=TranSectProfile(fn,ll,gr,varname)

%Linear Transect

addpath PlotProfile/
%zoom to region of inters


h = helpdlg('Zoom To Region Of interest. Klick ok when ready');

uiwait(h)
figure(1)

%place line
hwait = helpdlg('Draw Line to Create Transect');
uiwait(hwait)
clear h
h = imline(gca);


h2 = helpdlg('Continue ?');
uiwait(h2)


%get beginning and end
c=get(h,'Children');

xtrans=get(c(end),'Xdata');
ytrans=get(c(end),'Ydata');



delete(h);


hold on
plot(xtrans,ytrans,'k','linewidth',2)



%Find COntaining Polygon for ends else search next neighbours
Ptrans(1,:)=[xtrans(1) ytrans(1)];
Ptrans(2,:)=[xtrans(2) ytrans(2)];

parents=find_parent(ll,Ptrans);


cLat=(ytrans(1)+ytrans(2))/2;
cLon=(xtrans(1)+xtrans(2))/2;
[X, Y]=ll2xy(ll.hgrid.y,ll.hgrid.x,cLat,cLon);

grLocal=gr;%local projection grid
grLocal.hgrid.nodes(:,2)=X;
grLocal.hgrid.nodes(:,3)=Y;
grLocal.hgrid.x=X;
grLocal.hgrid.y=Y;

%Next neighbour for beining end end
if parents(1)==0
    [~, wo(1)]=min((ll.hgrid.x-xtrans(1)).^2+(ll.hgrid.y-ytrans(1)).^2);
    Ptrans(1,:)=[gr.hgrid.x(wo(1)) gr.hgrid.y(wo(1))];
    
    %translate to x,y coordinates of hgrid
    %P1=[gr.hgrid.x(wo(1))  gr.hgrid.y(wo(1))];
    P1=[X(wo(1))  Y(wo(1))];
    
else
    

    [P1(1), P1(2)]=ll2xy(Ptrans(1,2),Ptrans(1,1),cLat,cLon);
    
end

if parents(2)==0
    [~, wo(2)]=min((ll.hgrid.x-xtrans(2)).^2+(ll.hgrid.y-ytrans(2)).^2);
    %Ptrans(2,:)=[gr.hgrid.x(wo(2)) gr.hgrid.y(wo(2))];
    P2=[X(wo(2))  Y(wo(2))];
else
    [P2(1), P2(2)]=ll2xy(Ptrans(2,2),Ptrans(2,1),cLat,cLon);
end


%Create Local Projection centered around middle of transect
if P1(1)>P2(1)  %sort increasing from left to right
    P3=P1;
    P1=P2;
    P2=P3;
end


%enter spacing

dp=sqrt(sum((P2-P1).^2));

prompt = {'Enter grid point spacing (in m) for interpolated transect:'};
dlg_title = 'Input for peaks function';
num_lines = 1;
def = {'200',};
answer = inputdlg(prompt,dlg_title,num_lines,def);

ds=str2double(def{1});
np=round(dp/ds); %number of points
ds=dp/np;

dxx=(P2(1)-P1(1))/np;
dyy=(P2(2)-P1(2))/np;
xs=P1(1)+dxx*(0:np-1);
ys=P1(2)+dyy*(0:np-1);



Ps=[xs' ys'];
parents=find_parent(grLocal,Ps);%find containing elements


ivalid=parents~=0;
iquads(ivalid)=gr.hgrid.elem(parents(ivalid),2)==4;
itris(ivalid)=gr.hgrid.elem(parents(ivalid),2)==3;



%calculate distance weight - %triangles
% dxs=gr.hgrid.x(gr.hgrid.elem(parents(itris),3:5))-repmat(Ps((itris),1),1,3);
% dys=gr.hgrid.y(gr.hgrid.elem(parents(itris),3:5))-repmat(Ps((itris),2),1,3);
dxs=X(gr.hgrid.elem(parents(itris),3:5))-repmat(Ps((itris),1),1,3);
dys=Y(gr.hgrid.elem(parents(itris),3:5))-repmat(Ps((itris),2),1,3);




dR=sqrt(dxs.^2+dys.^2);
weights.tris=1./dR./repmat(sum(1./dR,2),1,3);
weights.tris(isnan(weights.tris))=1;

%calculate distance weight - %quadrangles
dxs=gr.hgrid.x(gr.hgrid.elem(parents(iquads),3:6))-repmat(Ps((iquads),1),1,4);
dys=gr.hgrid.y(gr.hgrid.elem(parents(iquads),3:6))-repmat(Ps((iquads),2),1,4);
dR=sqrt(dxs.^2+dys.^2);
weights.quad=1./dR./repmat(sum(1./dR,2),1,4);
weights.quad(isnan(weights.quad))=1;



%tangente
tangent=P2-P1;
nor=[tangent(2) -tangent(1)];
tangent=tangent/norm(tangent);
nor=nor/norm(nor);





hold on
Ar=diff(ylim)./diff(xlim); %aspect ratio
cx=mean(xtrans);
cy=mean(ytrans);
quiver(cx,cy,tangent(1)/100*Ar,tangent(2)*cos(cy*pi/180)/100,'k')    
quiver(cx,cy,nor(1)/100*Ar,nor(2)*cos(cy*pi/180)/100,'k')    


% 
% figure
% gr_plot(gr.hgrid,gr.hgrid.depth);
% axis([min(Ps(:,1)) max(Ps(:,1)) min(Ps(:,2)) max(Ps(:,2))])
% hold on
% plot(Ps(:,1),Ps(:,2),'k','linewidth',2)
% cx=mean(Ps(:,1));
% cy=mean(Ps(:,2));
% hold on
% quiver(cx,cy,tangent(1)*1000,tangent(2)*1000,'w','linewidth',2)    
% quiver(cx,cy,nor(1)*200,nor(2)*200,'w','linewidth',2)    



%perfekt match  1/0 -> inf   | inf/inf=nan  | nans have to get weight 1
interpVal=zeros(np,1);

%Interpol Distance weighted
interpVal(itris)=sum(gr.hgrid.depth(gr.hgrid.elem(parents(itris),3:5)).*weights.tris,2);
interpVal(iquads)=sum(gr.hgrid.depth(gr.hgrid.elem(parents(iquads),3:6)).*weights.quad,2);

bottom=interpVal;

s=(0:np-1)*ds;




tris=gr.hgrid.elem(parents(itris),3:5);
quads=gr.hgrid.elem(parents(iquads),3:6);

stack0=1;
stack=2;

stack=10;

fname=sprintf('%i_zcor.63',stack);
hzcor=sz_readHeader(fullfile(fn.outDir,fname));
nLevels=hzcor.vgrid.nLevels;


hvert=hzcor;
fname=sprintf('%i_vert.63',stack);
hvert.fname=fullfile(fn.outDir,fname);


fname2=sprintf('%i_%s',stack,varname);
hdata=sz_readHeader(fullfile(fn.outDir,fname2));



fname2=sprintf('%i_hvel.64',stack);
hhvel=sz_readHeader(fullfile(fn.outDir,fname2));



%load data;
zcor=sz_readTimeStep(hzcor,1);
zcor=reshape(zcor,hzcor.vgrid.nLevels,hzcor.hgrid.np);

if varname(end-1:end)=='64'
   keyboard 
end

data=sz_readTimeStep(hdata,1);
data=reshape(data,hzcor.vgrid.nLevels,hdata.hgrid.np);


vert=sz_readTimeStep(hvert,1);
vert=reshape(vert,hzcor.vgrid.nLevels,hzcor.hgrid.np);

hvel=sz_readTimeStep(hhvel,1);

u=reshape(hvel(:,1),nLevels,hzcor.hgrid.np);
v=reshape(hvel(:,2),nLevels,hzcor.hgrid.np);



dmesh=zeros(np,hzcor.vgrid.nLevels);
interpData=dmesh;
VV=dmesh;
UU=dmesh;
WW=dmesh;

for isigma=1:hzcor.vgrid.nLevels
    
    slayer=zcor(isigma,:);
    dmesh(itris,isigma)=sum(slayer(tris).*weights.tris,2);
    dmesh(iquads,isigma)=sum(slayer(quads).*weights.quad,2);
    
    %interpData
    dataLayer=data(isigma,:);
    dataLayer(dataLayer==-99)=nan;
    
    interpData(itris,isigma)=sum(dataLayer(tris).*weights.tris,2);
    interpData(iquads,isigma)=sum(dataLayer(quads).*weights.quad,2);
    
    %vertical velocity component
    dataLayer=vert(isigma,:);
    dataLayer(dataLayer==-9999)=nan;
    WW(itris,isigma)=sum(dataLayer(tris).*weights.tris,2);
    WW(iquads,isigma)=sum(dataLayer(quads).*weights.quad,2);
    
    
    % Transect Parelel end orthognal velocity
   dataLayer=u(isigma,:);
   dataLayer(dataLayer==-9999)=nan;
   UU(itris,isigma)=sum(dataLayer(tris).*weights.tris,2);
   UU(iquads,isigma)=sum(dataLayer(quads).*weights.quad,2);
   
   dataLayer=v(isigma,:);
   dataLayer(dataLayer==-9999)=nan;
   VV(itris,isigma)=sum(dataLayer(tris).*weights.tris,2);
   VV(iquads,isigma)=sum(dataLayer(quads).*weights.quad,2);

   
   %Rotate basis to unitvectors of transent tagent and normal
   nb=basis(tangent,nor,UU(:,isigma),VV(:,isigma));
    
   UU(:,isigma)=nb(1,:); %now Along channel Velocity
   VV(:,isigma)=nb(2,:);%now Across channel Velocity
     
    
end


dmesh(:,1)=-bottom;
dry=dmesh(:,2)==dmesh(:,end);
dmesh(dry,:)=repmat(bottom(dry),1,nLevels);



bottomMat=repmat(bottom,1,nLevels);
ifix=dmesh>bottomMat ;
dmesh(ifix)=bottomMat(ifix);




XX=repmat(s',1,nLevels);
YY=dmesh;

%UU=nan(size(XX));
%VV=UU;
%CC=dmesh;
%NN=UU;

%trasnect structure
transect.grLocal=grLocal;
transect.tangent=tangent;
transect.nor=nor;
transect.bottom=bottom;
transect.tris=tris;
transect.quads=quads;
transect.itris=itris;
transect.iquads=iquads;
transect.XX=XX;
transect.weights=weights;
transect.X=xs;
transect.Y=ys;
[transect.Lat,transect.Lon]=xy2ll(xs,ys,cLat,cLon);


figure
[A, Q]=profile_section(XX,YY,double(UU),double(WW)*5,interpData,VV);
colorbar
