function hovmuellerProfile(fn,ll,hovmue,dummy1,dummy2,Dim,nrvt,varname)

% function to visualizes of hovmueller diagrams
% of profile at point location
% - finds parent grid

fs=14; % fontsize used in plots

hold on
plot(hovmue.lon,hovmue.lat,'rx','linewidth',2)





%% parents - weights
parents=find_parent(ll,[hovmue.lon,hovmue.lat]);


elems=ll.hgrid.elem(parents,3:6);
x=hovmue.lon;
y=hovmue.lat;

x_x1=[ll.hgrid.x(elems(1)),x]';
x_x1=x_x1(:);
y_y1=[ll.hgrid.y(elems(1)),y]';
y_y1=y_y1(:);
x_x2=[ll.hgrid.x(elems(2)),x]';
x_x2=x_x2(:);
y_y2=[ll.hgrid.y(elems(2)),y]';
y_y2=y_y2(:);
x_x3=[ll.hgrid.x(elems(3)),x]';
x_x3=x_x3(:);
y_y3=[ll.hgrid.y(elems(3)),y]';
y_y3=y_y3(:);


d1=m_lldist(x_x1,y_y1);d1=d1(1:2:end);
d2=m_lldist(x_x2,y_y2);d2=d2(1:2:end);
d3=m_lldist(x_x3,y_y3);d3=d3(1:2:end);


d4=inf(size(d3));

for i=1:length(x)
    if ~isnan(elems(i,4))
        d4(i)=m_lldist([ll.hgrid.x(elems(i,4)),x(i)],...
            [ll.hgrid.y(elems(i,4)),y(i)]);
    end
end


%estimate node coordinate in gr3 dimensions
d_inv_tot=1./d1+1./d2+1./d3+1./d4;

w(:,1)=1./d1./d_inv_tot;
w(:,2)=1./d2./d_inv_tot;
w(:,3)=1./d3./d_inv_tot;
w(:,4)=1./d4./d_inv_tot;

ivalid=~isnan(elems);
depth=dot(w(ivalid),ll.hgrid.depth(elems(ivalid))');



%% Loda Data
dimname=sprintf('%s.%i',varname,Dim);
stack=hovmue.istack1;
hdata=sz_readHeader(fullfile(fn.outDir,sprintf('%i_%s',stack,dimname)));
if Dim==63
    hzcor=hdata;
    hzcor.fname=fullfile(fn.outDir,sprintf('%i_zcor.63',stack));
else
    hzcor=sz_readHeader(fullfile(fn.outDir,sprintf('%i_zcor.63',stack)));
end


ti=1:hdata.nSteps;
nt=length(ti);




Dinterp=zeros(nrvt,nt,'single');  %dummy second coordinate to  be able using interp2 for interpolation
weights.P1=repmat(w(1),nrvt,nt);
weights.P2=repmat(w(2),nrvt,nt);
weights.P3=repmat(w(3),nrvt,nt);
if ivalid(4)
    weights.P4=repmat(w(4),nrvt,nt);
end


ntot=(hovmue.istack2-hovmue.istack1+1)*nt;
Ztot=zeros(nrvt,ntot);
Dtot=Ztot;
timeinds=ti;

iup=3:nrvt;
idown=1:nrvt-2;

% for tris

if ~ivalid(4)
    
    for stack=hovmue.istack1:hovmue.istack2
    
    hzcor.fname=fullfile(fn.outDir,sprintf('%i_zcor.63',stack));
    hdata.fname=fullfile(fn.outDir,sprintf('%i_%s',stack,dimname));    
    
        
    data=sz_readTimeStep(hdata,ti);
    zcor=sz_readTimeStep(hzcor,ti);
        
    % parent element nodes sigmalevel timeseires
    %Zcor
    % Data
    inds=(elems(1)-1)*nrvt+(1:nrvt);
    Zparent1=zcor(inds,:);
    inds=(elems(2)-1)*nrvt+(1:nrvt);
    Zparent2=zcor(inds,:);
    inds=(elems(3)-1)*nrvt+(1:nrvt);
    Zparent3=zcor(inds,:);

    
    Z=Zparent1.*weights.P1+Zparent2.*weights.P2+Zparent3.*weights.P3;
    
    
    % Data
    inds=(elems(1)-1)*nrvt+(1:nrvt);
    Dparent1=data(inds,:);
    inds=(elems(2)-1)*nrvt+(1:nrvt);
    Dparent2=data(inds,:);
    inds=(elems(3)-1)*nrvt+(1:nrvt);
    Dparent3=data(inds,:);
    
    
    %liniear vertical interpolation to zlevels of profile location
    
    % surface and bottom are along sigma layer
    Dinterp([1 end],:)=Dparent1([1 end],:).*weights.P1([1 end],:)+...
        Dparent2([1 end],:).*weights.P2([1 end],:)+...
        Dparent3([1 end],:).*weights.P3([1 end],:);
    
    
    %linear vertical interpolation
    
    %Node 1
    intpData1=verticalInterp(Zparent1,Z,Dparent1);
    intpData2=verticalInterp(Zparent2,Z,Dparent2);
    intpData3=verticalInterp(Zparent2,Z,Dparent3);
    
    
    %horizontal inverse distance weighted interpolation on zlevel interpolated
    %data  | use iup only for size reasons weight matrixs is vertically equal
    Dinterp(2:end-1,:)=intpData1.*weights.P1(iup,:)+...
        intpData2.*weights.P2(iup,:)+...
        intpData3.*weights.P3(iup,:);
    
    
    Ztot(:,timeinds)=Z;
    Dtot(:,timeinds)=Dinterp;
    
    timeinds=timeinds+nt;
    
    end
    
    timemat=repmat(1:ntot,nrvt,1)';
    Ztot=Ztot';
    Dtot=Dtot';
    
    figure
    patchprofile(ntot,hdata.vgrid.nLevels,timemat(:),Ztot(:),Dtot(:));
    axis tight
    xlabel('time / h')
    ylabel('Depth / m')
    ch=colorbar;
    set(get(ch,'Ylabel'),'string',varname);
    
    
    
else quads
    
    keyboard
    % parent element nodes sigmalevel timeseires
    %Zcor
    % Data
    inds=(elems(1)-1)*nrvt+(1:nrvt);
    Zparent1=zcor(inds,:);
    inds=(elems(2)-1)*nrvt+(1:nrvt);
    Zparent2=zcor(inds,:);
    inds=(elems(3)-1)*nrvt+(1:nrvt);
    Zparent3=zcor(inds,:);
    inds=(elems(4)-1)*nrvt+(1:nrvt);
    Zparent4=zcor(inds,:);
    
    Z=Zparent1.*weights.P1+Zparent2.*weights.P2+Zparent3.*weights.P3+Zparent4.*weights.P4;
    
    
    
    % Data
    inds=(elems(1)-1)*nrvt+(1:nrvt);
    Dparent1=data(inds,:);
    inds=(elems(2)-1)*nrvt+(1:nrvt);
    Dparent2=data(inds,:);
    inds=(elems(3)-1)*nrvt+(1:nrvt);
    Dparent3=data(inds,:);
    inds=(elems(4)-1)*nrvt+(1:nrvt);
    Dparent4=data(inds,:);
    
    
    
    
    
    
    
end
% calc sigmalyers at profile location


%make Data available in Workspace
varname(varname=='-')='_';
eval(['workSpaceExport.hovmue.' varname  '=hovmue' ]);
%assignin('base','result',{localData,localData})
assignin('base','result',workSpaceExport)



    function intpData=verticalInterp(Zin,Zout,Din)
        
        Dup=Zin(iup,:)-Zout(2:end-1,:);
        Ddown=Zout(2:end-1,:)-Zin(idown,:);
        Dinv_tot=1./Dup+1./Ddown;
        
        
        % weights
        Wup=1./Dup./Dinv_tot;
        Wdown=1-Wup;
        
        intpData=Wup.*Din(iup,:)+Wdown.*Din(idown,:);
        
    end







end

