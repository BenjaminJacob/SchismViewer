
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

keyboard

%Find COntaining Polygon for ends else search next neighbours
Ptrans=[xtrans ytrans];
parents=find_parent(ll,Ptrans);

%Make Local Projection
cLat=mean(Ptrans(:,1));
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
ds=str2double(def{1});
%np=round(dp/ds); %number of points
%ds=dp/np;
L=sum(dp);
np=round(L/ds);
t=(0:np)*ds/L;
t(end)=min(t(end),1);
[pt,dudt,fofthandle] = interparc(t,P(:,1),P(:,2));

% 
% %sqrt(sum(diff(pt).^2,2));
% 
% 
% dxx=(P2(1)-P1(1))/np;
% dyy=(P2(2)-P1(2))/np;
% xs=P1(1)+dxx*(0:np-1);
% ys=P1(2)+dyy*(0:np-1);
% 

xs=pt(:,1);
ys=pt(:,2);

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
hovmue.x=xs;
hovmue.y=ys;
[hovmue.lon, hovmue.lat]=xy2ll(xs,ys,cLat,cLon);


interpVal(itris)=sum(gr.hgrid.depth(gr.hgrid.elem(parents(itris),3:5)).*weights.tris,2);
interpVal(iquads)=sum(gr.hgrid.depth(gr.hgrid.elem(parents(iquads),3:6)).*weights.quad,2);
hovmue.bottom=interpVal;
%plot(hovmue.bottom);


figure(1)
delete(ph1)
hold on
plot(hovmue.lon,hovmue.lat,'r');
 
%%

if XYZsource~=1
    
    %convert from lon lat to gr3
    elems=gr.hgrid.elem(:,3:6);
    
    
    
    x_x1=[ll.hgrid.x(elems(simplexid,1)),x]';
    x_x1=x_x1(:);
    y_y1=[ll.hgrid.y(elems(simplexid,1)),y]';
    y_y1=y_y1(:);
    x_x2=[ll.hgrid.x(elems(simplexid,2)),x]';
    x_x2=x_x2(:);
    y_y2=[ll.hgrid.y(elems(simplexid,2)),y]';
    y_y2=y_y2(:);
    x_x3=[ll.hgrid.x(elems(simplexid,3)),x]';
    x_x3=x_x3(:);
    y_y3=[ll.hgrid.y(elems(simplexid,3)),y]';
    y_y3=y_y3(:);
    
    
    
    d1=m_lldist(x_x1,y_y1);d1=d1(1:2:end);
    d2=m_lldist(x_x2,y_y2);d2=d2(1:2:end);
    d3=m_lldist(x_x3,y_y3);d3=d3(1:2:end);
    
    
    d4=inf(size(d3));
    
    for i=1:length(x)
        if ~isnan(elems(simplexid(i),4))
            d4(i)=m_lldist([ll.hgrid.x(elems(simplexid(i),4)),x(i)],...
                [ll.hgrid.y(elems(simplexid(i),4)),y(i)]);
            
        end
    end
    
    
    
       
    %estimate node coordinate in gr3 dimensions
    d_inv_tot=1./d1+1./d2+1./d3+1./d4;
    
    w(:,1)=1./d1./d_inv_tot;
    w(:,2)=1./d2./d_inv_tot;
    w(:,3)=1./d3./d_inv_tot;
    w(:,4)=1./d4./d_inv_tot;
    
    
    elems=elems(simplexid,:);
    elems(isnan(elems))=1;
    
    if size(elems,1)==1
        
        x=dot(w,gr.hgrid.x(elems)',2);
        y=dot(w,gr.hgrid.y(elems)',2);
        
    else
        
        x=dot(w,gr.hgrid.x(elems),2);
        y=dot(w,gr.hgrid.y(elems),2);
        
        
    end
    
    npoints=length(x)*length(zvals);
    XYZ=zeros(npoints,4);
    XYZ(:,1)=1:npoints;
    for i=1:length(zvals)
        XYZ((i-1)*length(x)+(1:length(x)),2:3)=[x y];
        XYZ((i-1)*length(x)+(1:length(x)),4)=zvals(i);
    end
    
    
    
end

