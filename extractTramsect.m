

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
plot(xtrans,ytrans,'r--')



%Find COntaining Polygon for ends else search next neighbours
Ptrans(1,:)=[xtrans(1) ytrans(1)];
Ptrans(2,:)=[xtrans(2) ytrans(2)];

parents=find_parent(gr,Ptrans);


%Next neighbour for beining end end
if parents(1)==0
    [was wo(1)]=min((ll.hgrid.x-xtrans(1)).^2+(ll.hgrid.y-ytrans(1)).^2);
    Ptrans(1,:)=[gr.hgrid.x(wo(1)) gr.hgrid.y(wo(1))];
    
    %translate to x,y coordinates of hgrid
    P1=[gr.hgrid.x(wo(1))  gr.hgrid.y(wo(1))];
    
else
    
    P1=Ptrans(1,:);
    
end

if parents(2)==0
    [was wo(2)]=min((ll.hgrid.x-xtrans(2)).^2+(ll.hgrid.y-ytrans(2)).^2);
    Ptrans(2,:)=[gr.hgrid.x(wo(2)) gr.hgrid.y(wo(2))];
    
    P2=[gr.hgrid.x(wo(2))  gr.hgrid.y(wo(2))];
else
    
    P2=Ptrans(2,:);
    
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


hold on
plot(xs,ys,'ro')

Ps=[xs' ys'];
parents=find_parent(gr,Ps);%find containing elements


iquads=gr.hgrid.elem(parents,2)==4;
itris=gr.hgrid.elem(parents,2)==3;



%calculate distance weight - %triangles
 dxs=gr.hgrid.x(gr.hgrid.elem(parents(itris),3:5))-repmat(Ps((itris),1),1,3);
 dys=gr.hgrid.y(gr.hgrid.elem(parents(itris),3:5))-repmat(Ps((itris),2),1,3);
 dR=sqrt(dxs.^2+dys.^2);
 weights.tris=1./dR./repmat(sum(1./dR,2),1,3);
weights.tris(isnan(weights.tris))=1; 
 
 %calculate distance weight - %quadrangles
 dxs=gr.hgrid.x(gr.hgrid.elem(parents(iquads),3:6))-repmat(Ps((iquads),1),1,4);
 dys=gr.hgrid.y(gr.hgrid.elem(parents(iquads),3:6))-repmat(Ps((iquads),2),1,4);
 dR=sqrt(dxs.^2+dys.^2);
 weights.quad=1./dR./repmat(sum(1./dR,2),1,4);
 weights.quad(isnan(weights.quad))=1; 
 
 
%perfekt match  1/0 -> inf   | inf/inf=nan  | nans have to get weight 1
interpVal=zeros(np,1);
 
%Interpol Distance weighted
interpVal(itris)=sum(gr.hgrid.depth(gr.hgrid.elem(parents(itris),3:5)).*weights.tris,2);
interpVal(iquads)=sum(gr.hgrid.depth(gr.hgrid.elem(parents(iquads),3:6)).*weights.quad,2);

bottom=interpVal;

s=(0:np-1)*ds;
figure
plot(s,interpVal)


figure
plot3(xs,ys,-interpVal)
axe=axis;


tris=gr.hgrid.elem(parents(itris),3:5);
quads=gr.hgrid.elem(parents(iquads),3:6);

stack0=1;
stack=2;

stack=10;

fname=sprintf('%i_zcor.63',stack);
hzcor=sz_readHeader(fullfile(fn.outDir,fname));


fname2=sprintf('%i_%s',stack,varname);
hdata=sz_readHeader(fullfile(fn.outDir,fname2));


%load data;
zcor=sz_readTimeStep(hzcor,1);
zcor=reshape(zcor,hzcor.vgrid.nLevels,hzcor.hgrid.np);

data=sz_readTimeStep(hdata,1);
data=reshape(data,hzcor.vgrid.nLevels,hdata.hgrid.np);


dmesh=zeros(np,hzcor.vgrid.nLevels);
interpData=dmesh;



for isigma=1:hzcor.vgrid.nLevels
    slayer=zcor(isigma,:);
    dataLayer=data(isigma,:);
    dataLayer(dataLayer==-99)=nan;
    dmesh(itris,isigma)=sum(slayer(tris).*weights.tris,2);
    dmesh(iquads,isigma)=sum(slayer(quads).*weights.quad,2);
    
    
    interpData(itris,isigma)=sum(dataLayer(tris).*weights.tris,2);
    interpData(iquads,isigma)=sum(dataLayer(quads).*weights.quad,2);
end
dmesh(:,1)=-bottom;

dry=dmesh(:,2)==dmesh(:,end);
dmesh(dry,:)=repmat(bottom(dry),1,hzcor.vgrid.nLevels);



bottomMat=repmat(bottom,1,nLevels);
ifix=dmesh>bottomMat ;
dmesh(ifix)=bottomMat(ifix);

figure
plot(dmesh)

nLevels=hzcor.vgrid.nLevels;

XX=repmat(s',1,nLevels);
YY=dmesh;

UU=nan(size(XX));
VV=UU;
CC=dmesh;
NN=UU;

figure
[A Q]=profile_section(XX,YY,UU,VV,interpData,NN)

