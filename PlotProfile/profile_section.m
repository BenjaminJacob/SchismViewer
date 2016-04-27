function [A, Q]=profile_section(XX,YY,UU,VV,CC,NN,varargin)

%nVerticalLayers=10;


%input
%      all m x n matrices
%     , where m is Dimension of section length
%       and n Dimension of section hight
%XX    each col has X value range (left to right in Plot)
%YY    each cols has depth profile according to the X column
%UU    Velocity Matrix UU(XX,YY) velocity in XX direction
%VV    Velocity in in YY direction
%CC    %is matric of same Dimension for color Plot e.g. section normal
%velocity salt ...
%NN section normal velocity

%Optional
%Dmin    %coloraxis min
%Dmax    %coloraxis max

% 
% m=size(XX,1);
% n=size(XX,2);%number of sigma layers



%Areas
dx=diff(XX(:,2:end));
hR=YY(2:end,2:end)-YY(2:end,1:end-1);
hL=YY(1:end-1,2:end)-YY(1:end-1,1:end-1);
hC=0.5*(hR+hL);
Ai=hC.*dx;

%Velocity SUmme Oben unten
NN2=NN(:,1:end-1)+NN(:,2:end);
Vi=(NN2(1:end-1,:)+NN2(2:end,:))/4;

%Flow
Qi=Ai.*Vi;
A=sum(Ai(:));
Q=sum(Qi(:));


%how='off';

%figure('visible','off')
% %figure(1)
% figure(1)
% clf
set(gca,'fontsize',14)


if max(XX(:,1)) >= 1000
    XX=XX/1000;
    labelx='transect length / km';
else
    labelx='transect length / m';
end
%ch=sms_patch_nodes(tri_pt,[],XX(:),YY(:),CC(:));

% keyboard
% ch=sms_patch_nodesstep(tri_pt,[],XX(:),YY(:),CC(:),step,Dmin,Dmax);

nxpoints=size(XX,1);
nsigma=size(XX,2);


%num1=find( ~isnan(YY(:,1)),1,'first');
%num2=find( ~isnan(YY(:,1)),1,'last');
%keyboard

% %kosmetic
% %for ii=1:num1-1
% 
% YY(num1-1,2:end)=YY(num1,2:end)
% %end
% 
% %for ii=num2+1:nxpoints
% YY(num2+1,2:end)=YY(num1,2:end)
% %end

%keyboard

patchprofile(nxpoints,nsigma,XX(:),YY(:),CC(:));

%%set(get(ch,'Ylabel'),'string','section normal vel. / ms^{-1}')
hold on
%plot(XX(1:m),YY,'w--')
xlabel(labelx)
ylabel('height / m')
plot(XX(:,1),YY(:,1),'k','linewidth',2)
xlim([0 XX(end)])
k=1;%0%00%00;
kv=1;





scale=0.2;
vscale=200;% factor for vertical component to horizontal

% interpolation
%keyboard
 Fu = TriScatteredInterp(XX(:),YY(:),UU(:));
 Fv = TriScatteredInterp(XX(:),YY(:),VV(:));

nx=ceil(size(XX,1)/4); 
ny=ceil(size(XX,2)/2); %pfeil alle 2 sigma schichten



x=linspace(XX(1),XX(end),nx);
y=linspace(min(YY(:)),max(YY(:)),ny);

[xm, ym]=meshgrid(x,y);

u=Fu(xm,ym);
v=Fv(xm,ym)*vscale;

%iq=
r=sqrt(u.^2+v.^2);
u=u./r;
v=v./r;

px=[XX(:,1); XX(end:-1:1,1)];
py=[YY(:,1); YY(end:-1:1,end)];

isin=inpolygon(xm,ym,px,py);
%vscale=200;
%quiver(xm,ym,u,v*vscale,scale,'k','linewidth',0)
nonan=~isnan(u) &(abs(u) >0 | abs(v) >0 ) &isin;
%qp=quiver(xm(nonan),ym(nonan),u(nonan),v(nonan),scale,'k','linewidth',0);
qp=quiver(xm(nonan),ym(nonan),u(nonan),v(nonan),scale,'k');

