% [A, Q]=iterateTransect(fn,transect,varname,stack0,stack1,'transpics')

%[A, Q]=CalcTransectFlux(fn2,transect,varname,stack0,stack1);
function  [A, Q]=CalcTransectFlux(fn,transect,varname,stack0,stack1,varargin)


%trasnect structure
%grLocal=transect.grLocal;
tangent=transect.tangent;
nor=transect.nor;
bottom=transect.bottom;
tris=transect.tris;
quads=transect.quads;
itris=transect.itris;
iquads=transect.iquads;
XX=transect.XX;
weights=transect.weights;

stack=stack0;

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

np=size(transect.XX,1);
YY=zeros(np,hzcor.vgrid.nLevels);
interpData=YY;
VV=YY;
UU=YY;

nt=(stack1-stack0+1)*hdata.nSteps;
Q=zeros(nt,1);
A=Q;

ttot=0;

startdate=datenum(2011,7,1);


dx=diff(XX(:,2:end));

for stack=stack0:stack1
    
    fname=sprintf('%i_zcor.63',stack);
    hzcor.fname=fullfile(fn.outDir,fname);
    
    
    
    hvert=hzcor;
    fname=sprintf('%i_vert.63',stack);
    hvert.fname=fullfile(fn.outDir,fname);
    
    
    fname2=sprintf('%i_%s',stack,varname);
    hdata.fname=fullfile(fn.outDir,fname2);
    
    fname2=sprintf('%i_hvel.64',stack);
    hhvel.fname=fullfile(fn.outDir,fname2);
    
    
    for ti=1:24
        
        ttot=ttot+1;
        
        %load data;
        zcor=sz_readTimeStep(hzcor,ti);
        zcor=reshape(zcor,hzcor.vgrid.nLevels,hzcor.hgrid.np);
        
        if varname(end-1:end)=='64'
            keyboard
        end
        
        data=sz_readTimeStep(hdata,ti);
        data=reshape(data,hzcor.vgrid.nLevels,hdata.hgrid.np);
        
        
        vert=sz_readTimeStep(hvert,ti);
        vert=reshape(vert,hzcor.vgrid.nLevels,hzcor.hgrid.np);
        
        hvel=sz_readTimeStep(hhvel,ti);
        
        u=reshape(hvel(:,1),nLevels,hzcor.hgrid.np);
        v=reshape(hvel(:,2),nLevels,hzcor.hgrid.np);
        
        
        
        for isigma=1:hzcor.vgrid.nLevels
            
            slayer=zcor(isigma,:);
            YY(itris,isigma)=sum(slayer(tris).*weights.tris,2);
            YY(iquads,isigma)=sum(slayer(quads).*weights.quad,2);
            
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
        
        
        YY(:,1)=-bottom;
        dry=YY(:,2)==YY(:,end);
        YY(dry,:)=repmat(bottom(dry),1,nLevels);
        
        
        
        bottomMat=repmat(bottom,1,nLevels);
        ifix=YY>bottomMat ;
        YY(ifix)=bottomMat(ifix);
        
        
        
        
%       XX=repmat(s',1,nLevels);
%        YY=YY;
        
        %UU=nan(size(XX));
        %VV=UU;
        %CC=YY;
        %NN=UU;

        
        


%Areas
hR=YY(2:end,2:end)-YY(2:end,1:end-1);
hL=YY(1:end-1,2:end)-YY(1:end-1,1:end-1);
hC=0.5*(hR+hL);
Ai=hC.*dx;

%Velocity SUmme Oben unten
NN2=VV(:,1:end-1)+VV(:,2:end);
Vi=(NN2(1:end-1,:)+NN2(2:end,:))/4;

%Flow
Qi=Ai.*Vi;
A(ttot)=sum(Ai(:));
Q(ttot)=sum(Qi(:));


        
%         [A(ttot), Q(ttot)]=profile_section(XX,YY,double(UU),double(WW),VV,VV);
%         datum=startdate+ttot*hzcor.dt/86400;
%         title(datestr(datum))
%         ch=colorbar;
%         set(get(ch,'Ylabel'),'string','U_{cross} / ms^{-1}','fontsize',15)
%         %caxis([-1 1]*2.5)
%         caxis([-1 1]*1)
%         colormap redblue
%         if nargin==6
%         picname=sprintf('%03i_transect',ttot);
%         print('-dpng',fullfile(outdir,picname),'-r 600')
%         end
    end
    
end


end