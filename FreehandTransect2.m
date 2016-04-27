
function hovmue=FreehandTransect(fn,gr,ll)
%% create transect line

%Linear Transect

addpath PlotProfile/
%zoom to region of inters


h = helpdlg('Zoom to region Of interest. Klick ok when ready');

uiwait(h)
figure(1)

%place line
hwait = helpdlg('Click points of the line to create transect. Finish with double click');
uiwait(hwait)
clear h
h = impoly(gca);


h2 = helpdlg('Continue ?');
uiwait(h2)


%get beginning and end
c=get(h,'Children');


xtrans=get(c(end),'Xdata');
ytrans=get(c(end),'Ydata');



delete(h);

hold on
ph1=plot(xtrans,ytrans,'k.-','linewidth',2);


%Find COntaining Polygon for ends else search next neighbours
Ptrans=[xtrans ytrans];
parents=find_parent(ll,Ptrans);

%Make Local Projection
cLat=mean(Ptrans(:,2));
cLon=mean(Ptrans(:,1));
[X, Y]=ll2xy(ll.hgrid.y,ll.hgrid.x,cLat,cLon);

grLocal=gr;%local projection grid
grLocal.hgrid.nodes(:,2)=X;
grLocal.hgrid.nodes(:,3)=Y;
grLocal.hgrid.x=X;
grLocal.hgrid.y=Y;

%Next neighbour for beining end end

P=zeros(size(Ptrans));

for i=1:size(parents,1)
    
    
    if parents(i)==0
        [~, wo(i)]=min((ll.hgrid.x-xtrans(i)).^2+(ll.hgrid.y-ytrans(i)).^2);
        Ptrans(i,:)=[gr.hgrid.x(wo(i)) gr.hgrid.y(wo(i))];
        
        %translate to x,y coordinates of hgrid
        %P1=[gr.hgrid.x(wo(1))  gr.hgrid.y(wo(1))];
        P(i,:)=[X(wo(i))  Y(wo(i))];
        
    else
        
        [P(i,1), P(i,2)]=ll2xy(Ptrans(i,2),Ptrans(i,1),cLat,cLon);
    end
    
end


dp=sqrt(sum(diff(P).^2,2));
prompt = {'Enter grid point spacing (in m) for interpolated transect:'};
dlg_title = 'Input for peaks function';
num_lines = 1;
def = {'200',};
answer = inputdlg(prompt,dlg_title,num_lines,def);
ds=str2double(answer{1});

L=sum(dp);
np=round(L/ds);
t=(0:np)*ds/L;
t(end)=min(t(end),1);
[pt,~,~] = interparc(t,P(:,1),P(:,2));


%sqrt(sum(diff(pt).^2,2));


% dxx=(P2(1)-P1(1))/np;
% dyy=(P2(2)-P1(2))/np;
% xs=P1(1)+dxx*(0:np-1);
% ys=P1(2)+dyy*(0:np-1);


xs=pt(:,1);
ys=pt(:,2);


%Find triangular Weights and distance to parent elements
Ps=[xs ys];
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


hovmue.weights=weights;
hovmue.grLocal=grLocal;
hovmue.tris=gr.hgrid.elem(parents(itris),3:5);
hovmue.quads=gr.hgrid.elem(parents(iquads),3:6);

% converte coordinates on local projection back to lon lat
hovmue.xLocal=xs;
hovmue.yLocal=ys;
hovmue.ds=[0; sqrt(diff(xs).^2+diff(ys).^2)];
[hovmue.lat, hovmue.lon]=xy2ll(xs,ys,cLat,cLon);


interpVal(itris)=sum(gr.hgrid.depth(gr.hgrid.elem(parents(itris),3:5)).*weights.tris,2);
interpVal(iquads)=sum(gr.hgrid.depth(gr.hgrid.elem(parents(iquads),3:6)).*weights.quad,2);
hovmue.bottom=interpVal;
%plot(hovmue.bottom);


figure(1)
delete(ph1)
hold on
plot(hovmue.lon,hovmue.lat,'r');

%%


