function [openbd landbd ]=plotBound2(grid)

%n_nodes - nr of nodes of open boundaries
%openbd - openboundaries list

%plot boundary


iph=0; %plothandle for boundaries

disp('loading open and closed boundaries')
fid=fopen(grid.fname,'r');


textline=fgets(fid);
pattern=' Number of open boundaries';

while isempty(strfind(textline,pattern))
    textline=fgets(fid);
end

n_opbd=textscan(textline,'%d ');
n_opbd=n_opbd{1};

 textline=fgets(fid);
n_opbd_nodes=textscan(textline,'%d ');
n_opbd_nodes=n_opbd_nodes{1};
 
hold on
%openbd=zeros(n_opbd_nodes,1);
n_nodes=zeros(1,n_opbd);
k=1;
crits=[];%bnd nodes which may node be deep enouqh
for i=1:n_opbd
    
    
    textline=fgets(fid);%Number of nodes for open bundary i
    temp=textscan(textline,'%d ');
    n_nodes(i)=temp{1};
    
   
    %pattern= sprintf(' = Number of nodes for open boundary',bound_segs(i));

    
    
    temp=textscan(fid,'%d \n',n_nodes(i));%read boundaries segments (nr of node specified by bound segs(i)
    %openbd(k:k-1+n_nodes(i))=temp{1};
    openbd{i}=temp{1};
    %k=k+n_nodes(i);

    
    iph=iph+1;
    
end



%% land boundaries

textline=fgets(fid);

n_ldbd=textscan(textline,'%d '); %number of land boundaries
n_ldbd=n_ldbd{1};

 textline=fgets(fid);
n_ldbd_nodes=textscan(textline,'%d ');
n_ldbd_nodes=n_ldbd_nodes{1};
 
hold on
%landbd=zeros(n_ldbd_nodes,1);
n_lnodes=zeros(1,n_ldbd);
k=1;
for i=1:n_ldbd
    
    
    
    textline=fgets(fid);%Number of nodes for open bundary i
    temp=textscan(textline,'%d ');
    n_lnodes(i)=temp{1}(1);
    
   
    %pattern= sprintf(' = Number of nodes for open boundary',bound_segs(i));

    
    
    temp=textscan(fid,'%d \n',n_lnodes(i));%read boundaries segments (nr of node specified by bound segs(i)
    %landbd(k:k-1+n_lnodes(i))=temp{1};
    landbd{i}=temp{1};
    k=k+n_lnodes(i);

    
   iph=iph+1;
   %phbound(iph)=plot(grid.x(temp{1}),grid.y(temp{1}),'r');
    
end



