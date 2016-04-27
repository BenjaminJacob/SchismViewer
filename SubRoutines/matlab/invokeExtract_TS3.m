

function [XYZ ics ifs istack1 istack2]=invokeExtract_TS(fn,gr,ll,fill_in,basedir,binary_file,bp_fileName,bp_file,days)

%% get XYZ(points)

ok=0;
while ~ok
    
    
    XYZsource=menu('XYZ points','1) read in from .gr3 file','2) enter vectors of x,y and z coordinates','3) choose xy graphically and enter z' );
    
    ok=1;
    switch XYZsource
        
        case 1
            
            [file folder]=uigetfile('*.gr3','open gr3 file');
            gr_BP.hgrid=gr_readHGrid(fullfile(folder,file));
            npoints=gr_BP.hgrid.np;
            XYZ=gr_BP.hgrid.nodes;
            
            
            
            
        case 2
            prompt = {'Enter (GEO) LON coordinate:','Enter LAT coordinate:','Enter Z coordinate: '};
            dlg_title = 'Select position';
            num_lines = 1;
            def = {'0','0','0'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            x=str2double(answer(1));
            y=str2double(answer(2));
            
            %x=sscanf(answer{1},'%f');
            %y=sscanf(answer{2},'%f');
            zvals=sscanf(answer{3},'%f');
            
            
            %hold on
            plot(x,y,'k--+')
            text(x,y,num2str((1:length(x))','P%i'));
            
            
            
            %simplexid=intriangle(ll.hgrid,[x y]);
            simplexid=find_parent(ll.hgrid,[x y]);
            if sum(isnan(simplexid))
                disp('point(s) not in domian please reselect')
                ok=0;
                continue
            end
            
            
        case 3
            
            
            
            hhelp=helpdlg('In main window Click Tools to zoom to region. After pressing ok elect Points by clicking. Press return to finish selection');
            waitfor(hhelp) % continue after ok
            figure(1)
            [x, y]=ginput();
            %            simplexid=intriangle(ll.hgrid,[x y]);

            simplexid=find_parent(ll,[x y]);
            if isnan(simplexid)
                disp('point not in domian please reselect')
                ok=0;
                %    continue
            end
            
            hold on
            plot(x,y,'k--+')
            text(x,y,num2str((1:length(x))','P%i'));
            
            
            
            %inputdialog
            prompt = {'enter (space sepearated) z value(s) for the locations. Z values applied to all xy pairs.'};
            dlg_title = 'Input Z coordinates';
            num_lines = 1;
            def = {'0'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            
            zvals=sscanf(answer{1},'%f');
            %zvals=input('enter z values for the locations. All Z values will be applied to all xy pairs. Enter multiple Values seperated by space ','s');
            %zvals=sscanf(zvals,'%f');
            
    end
    
    
end



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


%hold on
%plot(XYZ(:,2),XYZ(:,3),'k--+')
%text(XYZ(:,2),XYZ(:,3),num2str((1:size(XYZ,1))','P%i'))



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


%% create extract optionsFile for fortran routine


% 
% fill_in=-9999;
% bp_file='station.bp';
% npoints=size(XYZ,1);
% 
% %% extract options
% fid=fopen('SubRoutines/extract_options','w');
% fprintf(fid,'fill_in=%i       !Fill in value used for abnormal case (e.g., below bottom, dry etc)\n',fill_in);
% fprintf(fid,'rundir="%s"   !directory where grid files are stored (e.ge one level above outputs)\n',fn.runDir);
% fprintf(fid,'basedir="%s"   !basedir - directory where binary and build point files are stored\n',fn.outDir);
% fprintf(fid,'binary_file=%s !binary_file - elev.61, hvel.64 etc; it should be inside basedir\n',binary_file);
% fprintf(fid,'istack1=%i !     istack[1,2] - start and end stack # for binary outputs;\n',stack1);
% fprintf(fid,'istack2=%i\n',stack2);
% fprintf(fid,'%s=%s !  %s - build point file (see station.ts.sample / station.xyzt.sample)\n','bp_file',bp_file,'bp_file');
% fprintf(fid,'%s=%i  ! (ics=1: Cartesian coordinates; =2: lat/lon.\n','ics',ics);
% fprintf(fid,'%s=%i !ifs=0: input z are z coordinates;ifs =1: input z are relative to free surface; and in this case, make sure all z>=0\n','ifs',ifs);
% fclose(fid)
% 
% 
% 
% 
% fid=fopen('SubRoutines/fortran/extract_options','w');
% fprintf(fid,'fill_in=%i       !Fill in value used for abnormal case (e.g., below bottom, dry etc)\n',fill_in)
% fprintf(fid,'rundir="%s"   !directory where grid files are stored (e.ge one level above outputs)\n',fn.runDir)
% fprintf(fid,'basedir="%s"   !basedir - directory where binary and build point files are stored\n',basedir)
% fprintf(fid,'binary_file=%s !binary_file - elev.61, hvel.64 etc; it should be inside basedir\n',binary_file)
% fprintf(fid,'istack1=%i !     istack[1,2] - start and end stack # for binary outputs;\n',istack1)
% fprintf(fid,'istack2=%i\n',istack2)
% fprintf(fid,'%s=%s !  %s - build point file (see station.ts.sample / station.xyzt.sample)\n',bp_fileName,bp_file,bp_fileName);
% fprintf(fid,'%s=%i  ! (ics=1: Cartesian coordinates; =2: lat/lon.\n','ics',ics);
% fprintf(fid,'%s=%i !ifs=0: input z are z coordinates;ifs =1: input z are relative to free surface; and in this case, make sure all z>=0\n','ifs',ifs);
% 
% 
% 
% %% create station.ts input file for fortran routine
% 
% 
% itransect=0; %0 time seires;  1: transect
% 
% %ics=2; %1: cartesian 2: lat lon
% %ics=input('are the coordinates 1) cartesian  or 2) in lat lon ? Enter 1 or 2: ');
% 
% %ifs=input('are input z coordinates 0) z cordinates or 1) relative to free surface; and in this case, make sure all z>=0 ? Enter 0 or 1: ');
% %ifs=0; %0: z are z coordinates; 1: z relative to free surface (nedd alll z>0)
% %folder='./outputs/';
% %system(sprintf('ls  %s', fnA.outDir))
% fid=fopen(fullfile('SubRoutines/fortran',bp_file),'w');
% fprintf(fid,'%i %i %i  !itransect,ics,ifs (itransect=0: time series; =1:transect. ics=1: Cartesian coordinates; =2: lat/lon. ifs=0: input z are z coordinates; =1: input z are relative to free surface; and in this case, make sure all z>=0 :\n',itransect,ics,ifs)
% fprintf(fid,'%i !total # of points\n',npoints);
% fprintf(fid,'%4i %14.6f %14.6f %4.4f !point #, x, y, z\n',XYZ(1,:));
% if npoints>1
%     fprintf(fid,'%4i %14.6f %14.6f %4.4f\n',XYZ(2:end,:)');
% end
% fclose(fid);
% 
% 
% % Construct a questdlg with three options
% 

end
