

function hovmue=ExtractHovmueller(fn,gr,ll,varname,days)

%% get XYZ(points)



%hovmue=FreehandTransect(fn,gr,ll);

%
%
%
%
% %keyboard
% if nargin>4 && extractType==2
%     varnames={varnames 'zcor.63'};
% else
%     varnames={varnames};
% end



     if varname(end)<'3'
        extractType=menu('extract options','Create Horizonzal Transect');
     else
         extractType=menu('extract options','Create Horizonzal Transect','Create Vertical Profile');
     end



if extractType==1;
    hovmue=FreehandTransect2(fn,gr,ll);
    hovmue.type='transect';
else
    
     
     hovmue.type='profile';
    %place line
    hwait = helpdlg('Click points to define coordinates of Profile');
    uiwait(hwait)
    clear h
    [hovmue.lon, hovmue.lat]=ginput(1);
    
end
ok=0;
%
%
% prompt = {'enter stack 1: ','enter stack 2: '};
% dlg_title = 'Input Coordinate time intervall';
% num_lines = 1;
% def = {num2str(days.min),num2str(days.max)};
% answer = inputdlg(prompt,dlg_title,num_lines,def);
% istack1=str2double(answer{1});
% istack2=str2double(answer{2});




%
% nSteps=h.nSteps;
% nt=(istack2-istack1+1)*nSteps;
% nx=length(hovmue.bottom);
% hovmue.data=zeros(nt,nx);
%
% hovmue.data
%
%
%
%
% tridata=zeros(nx,3,nSteps);
% triweights=repmat(hovmue.weights.tris,1,1,nSteps);
%
%
% inds=1:nSteps;


% % %% stack 1   - weights method - Problem no vertical interpolations so far
%  stack=istack1;
%  h=sz_readHeader(fullfile(fn.outDir,sprintf('%i_%s',stack,varname)));
%
% ti=1:h.nSteps;
% subinds1=1:3:3*nx;
% subinds2=subinds1+1;
% subinds3=subinds2+1;
%
%
% data=sz_readTimeStep(h,ti);
% temp=data(hovmue.tris',1,:);
% tridata(:,1,:)=temp(subinds1,1,:);
% tridata(:,2,:)=temp(subinds2,1,:);
% tridata(:,3,:)=temp(subinds3,1,:);
%
% hovmue.data(inds,:)=squeeze(dot(tridata,triweights,2))';
%
% inds=inds+nSteps;
%
% %% stack2 - end
% for stack=istack1+1:istack2
%
%     h.fname=fullfile(fn.outDir,sprintf('%i_%s',stack,varname));
%
%     data=sz_readTimeStep(h,ti);
%     temp=data(hovmue.tris',1,:);
%     tridata(:,1,:)=temp(subinds1,1,:);
%     tridata(:,2,:)=temp(subinds2,1,:);
%     tridata(:,3,:)=temp(subinds3,1,:);
%
%     hovmue.data(inds,:)=squeeze(dot(tridata,triweights,2))';
%
%
%     inds=inds+nSteps;
% end
%
% figure
% imagesc(hovmue.data)
% xlabel('Position')
% ylabel('time')

%% Extract using fortran scripts



prompt = {'are input z coordinates 0) z cordinates or 1) relative to free surface; and in this case, make sure all z>=0 ?',...
    'enter stack 1: ','enter stack 2: '};
dlg_title = 'Input Coordinate type & time intervall';
num_lines = 1;
def = {'0',num2str(days.min),num2str(days.max)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
ics=1+gr.hgrid.x(1)==ll.hgrid.x(1);
ifs=str2double(answer{1});
istack1=str2double(answer{2});
istack2=str2double(answer{3});
%istack1=input('enter stack 1: ');
%istack2=input('enter stack 2: ');


if extractType==1;
    %inputdialog
    prompt = {'enter z value for Hovmueller diagram'};
    dlg_title = 'Input Z coordinates';
    num_lines = 1;
    def = {'0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    zvals=sscanf(answer{1}(1),'%f');
    
    
    x=hovmue.lon;
    y=hovmue.lat;
    
    hovmue.z=zvals;
    
    
    if gr.hgrid.x(1)~=ll.hgrid.x(1)
        
        %convert from lon lat to gr3
        elems=gr.hgrid.elem(:,3:6);
        
        
        
        
        x_x1=[ll.hgrid.x(elems(parents,1)),x]';
        x_x1=x_x1(:);
        y_y1=[ll.hgrid.y(elems(parents,1)),y]';
        y_y1=y_y1(:);
        x_x2=[ll.hgrid.x(elems(parents,2)),x]';
        x_x2=x_x2(:);
        y_y2=[ll.hgrid.y(elems(parents,2)),y]';
        y_y2=y_y2(:);
        x_x3=[ll.hgrid.x(elems(parents,3)),x]';
        x_x3=x_x3(:);
        y_y3=[ll.hgrid.y(elems(parents,3)),y]';
        y_y3=y_y3(:);
        
        
        
        d1=m_lldist(x_x1,y_y1);d1=d1(1:2:end);
        d2=m_lldist(x_x2,y_y2);d2=d2(1:2:end);
        d3=m_lldist(x_x3,y_y3);d3=d3(1:2:end);
        
        
        
        d4=inf(size(d3));
        
        for i=1:length(x)
            if ~isnan(elems(parents(i),4))
                d4(i)=m_lldist([ll.hgrid.x(elems(parents(i),4)),x(i)],...
                    [ll.hgrid.y(elems(parents(i),4)),y(i)]);
                
            end
        end
        
        
        
        
        %estimate node coordinate in gr3 dimensions
        d_inv_tot=1./d1+1./d2+1./d3+1./d4;
        
        w(:,1)=1./d1./d_inv_tot;
        w(:,2)=1./d2./d_inv_tot;
        w(:,3)=1./d3./d_inv_tot;
        w(:,4)=1./d4./d_inv_tot;
        
        
        elems=elems(parents,:);
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
        
        
    else % All in Lon Lat coordinates
        
        
        
        
        %hovmue.lon
        npoints=length(x)*length(zvals);
        
        
        XYZ=zeros(npoints,4);
        XYZ(:,1)=1:npoints;
        for i=1:length(zvals)
            XYZ((i-1)*length(x)+(1:length(x)),2:3)=[x y];
            XYZ((i-1)*length(x)+(1:length(x)),4)=zvals(i);
        end
        
        
        
        
        
        
    end
    
    
    
    choiceSave = questdlg('Would you like to save the coordinates used for extraction', ...
        'Save coordinates', ...
        'Yes','No','No');
    
    if strcmp(choiceSave,'Yes')
        [FileName,PathName] = uiputfile('*.gr3');
        fname=(fullfile(PathName,FileName));
        fid=fopen(fname,'w');
        fprintf(fid,'%s\n',fname);
        fprintf(fid,'%i\n',npoints);
        fprintf(fid,'%u %6.6f %6.6f %g \n',XYZ');
    end
    
    
    hovmue.XYZ=XYZ;
    hovmue.ics=ics;
    hovmue.ifs=ifs;
    
    
end
hovmue.istack1=istack1;
    hovmue.istack2=istack2;
end
