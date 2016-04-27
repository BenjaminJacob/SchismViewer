function SCHISM_Data_Inspector5
% SIMPLE_GUI2 Select a data set from the pop-up menu, then
% click one of the plot-type push buttons. Clicking the button
% plots the selected data in the axes.
close all


%% add pathes

addpath(genpath('thirdPartyToolbox')); % add toolboxes used like m-elio for plotting routines
addpath(genpath('thirdData')); %add Bathymetry Datta and coastlines
gshhs_file='gshhg-bin-2.3.4/gshhs_f.b';
addpath(genpath('SubRoutines')); %add subroutines



% Plot Boundaries from shapefile from
% http://www.naturalearthdata.com/
%useShapefile=1;
%shapefile=[pwd '/ne_50m_land/ne_50m_land.shp' ];
% if useShapefile
%     landplt=geoshow(shapefile);
%     set(findall(landplt,'type','patch'),'FaceColor',[.7 .7 .7]);
%     axis(axe)
% end


%addpath(genpath(pwd))
%% Initilaize with file pathes

%add pathes
addpath(genpath('SubRoutines'));



% open run Setupfile
fn.runDir = uigetdir([],'Select Selfe run directory - where hgrid.gr3 ');

% open combined file
fn.outDir = uigetdir(fn.runDir,'Select directory of combined binaries - eg 1_elev.61 ');

fn2=[];


%creat link copy in fortran source directory using lndir
%fn.outDir
if exist('outputs_linked','dir')
    !rm -r outputs_linked
end
!mkdir outputs_linked
system(sprintf('lndir %s outputs_linked',fn.outDir))


% Global Variables
%gridName='hgrid.gr3';
gridName='hgrid.gr3';
disp('loading hgrid.gr3');gr.hgrid=gr_readHGrid(fullfile(fn.runDir,gridName));
gr2=[];
%hgird ll
gridName='hgrid.ll';
disp('loading hgrid.ll');ll.hgrid=gr_readHGrid(fullfile(fn.runDir,gridName));
assignin('base','gr',gr)
assignin('base','ll',ll)


dummy1=-99;
dummy2=-9999;
alternate=0; %alternating switch to plot depth and depth difference altenating
inregnodes=1:gr.hgrid.np;


% count stacks
elevFiles=dir([fn.outDir '/*elev.61' ]);
stackList=cell(length(dir([fn.outDir '/*elev.61' ]))+1,1);
stackList{1}='stack';
for j=1:length(stackList)-1;
    stacknr=elevFiles(j).name(1:find(elevFiles(j).name=='_')-1);
    stackList{j+1}=stacknr;
end

for i=2:length(stackList)
    stackList{i}=str2double(stackList{i});
end
temp=(sort(cell2mat(stackList(2:end))));
for i=2:length(stackList)
    stackList{i}=num2str(temp(i-1));
end

days.min=stackList{2};
days.max=stackList{end};


%varlist Different variables in folder
varlist=dir([fn.outDir '/' stackList{2} '_*.6*' ]);
varNames=cell(length(varlist)+1,1);
varDims=zeros(length(varlist)+1,1);

varNames{1}='depth';
varDims(1)=0;

for j=1:length(varlist);
    iunScore=find(varlist(j).name=='_');
    idot=find(varlist(j).name=='.');
    varNames{j+1}= varlist(j).name(iunScore+1:idot-1);
    varDims(j+1)=str2double(varlist(j).name(idot+1:end));
end



handles.tooltip.mesh=0;


timesteplist={};
timesteplist{1}='timestep';
difference_plot=0;
ts1=[];
ts2=[];
dts=[];

transect=struct;
%placeholders
video=[];
hovmue=[];
%video.Patchdata
%video.stack0
%video.stack1
%% Create Gui
posX=360;
posY=500;

%get(0,'MonitorPositions')
scrsz = get(0,'ScreenSize');

LY=scrsz(end)*0.5;
LX=round(1.41*LY);


xDrop=2;
%eLEMENT WIDTH
%hx=100;
hy=25;
%hxdrop=80;
hxdrop=([250 60 80 120]);
yDrop=LY-22-hy;
%
%  Create and then hide the GUI as it is being constructed.
f = figure('Visible','off','Position',[posX,posY,LX,LY],'KeyPressFcn',@keypress);
% Assign the GUI a name to appear in the window title.
set(f,'Name','SCHISM View Tool')
% Move the GUI to the center of the screen.
movegui(f,'center')

% Make the GUI visible.
set(f,'Visible','on');
% Construct the components.

%hpoint    = uicontrol('Style','pushbutton',...
%   'String','Select Point','Position',[LX-70*3,LY-hy,70*2,hy],...
%  'Callback',{@SelectPoint_Callback});
hplot    = uicontrol('Style','pushbutton',...
    'String','plot','Position',[LX-70,LY-hy,70,hy],...
    'Callback',{@Plotbutton_Callback});

hextract    = uicontrol('Style','pushbutton',...
    'String','Extract ...','Position',[LX-70*3,LY-hy,70*2,hy],...
    'Callback',{@ExtractButton_Callback},'ToolTipString','Extract either time series, slabs, or xyzt ');


hcreate    = uicontrol('Style','pushbutton',...
    'String','Create ...','Position',[LX-70*3,LY-2*hy,70*2,hy],...
    'Callback',{@CreateButton_Callback},'ToolTipString','Extract either time series, slabs, or xyzt ');


hvideo = uicontrol('Style','pushbutton',...
    'String','video','Position',[LX-70,LY-2*hy,70,hy],...
    'Callback',{@VideoButton_Callback});

hexport    = uicontrol('Style','pushbutton',...
    'String','export','Position',[LX-70,LY-3*hy,70,hy],...
    'Callback',{@exportButton_Callback});





% htext  = uicontrol('Style','text','String','Select Data',...
%     'Position',[325,90,60,15]);

hlabel = uicontrol('Style','text',...
    'String','select: stack and timestep to load data',...
    'Position',[1,LY-hy-1,hxdrop(3)*3.5,hy]);

hsetup    = uicontrol('Style','pushbutton',...
    'String','Select Setup','Position',[hxdrop(3)*3.5+2,LY-hy,70*2,hy],...
    'Callback',{@SetupButton_Callback});



hpopup = uicontrol('Style','popupmenu',...
    'String',  varNames,...
    'Position',[xDrop,yDrop,hxdrop(1),hy],...
    'Callback',{@popup_menu_Callback});

hcheck = uicontrol('Style','checkbox',...
    'String',  'Difference between Setups (substract 2. from 1. selected)',...
    'Position',[xDrop,yDrop-hy*3/4,hxdrop(3)*5,hy*3/4],...
    'Callback',{@checkbox_Callback});



hpopupStack = uicontrol('Style','popupmenu',...
    'String',stackList,...
    'Position',[xDrop+hxdrop(1),yDrop,hxdrop(2),hy],...
    'Callback',{@popup_menu_Callback_Stack});




hpopupTime = uicontrol('Style','popupmenu',...
    'String',timesteplist,...
    'Position',[xDrop+sum(hxdrop(1:2)),yDrop,hxdrop(3),hy],...
    'Callback',{@popup_menu_Callback_Time});


%read number of layers from vgrid
found=1;
if exist(fullfile(fn.runDir, 'vgrid.in'),'file')
    fid=fopen(fullfile(fn.runDir, 'vgrid.in'),'r');
    line=fgets(fid);
    while isempty(strfind(line,'nvrt'))
        line=fgets(fid); %layer info in second linei
        if feof(fid)
            found=0;
            break
        end
    end
    
    
    if found
        nrvt=textscan(line,'%d');
        nrvt=nrvt{1}(1); %nr of vertical layers
        fclose(fid);
        
        
    else
        disp('vgrid.in does not excist')
        nrvt=input('enter # of vertical layers ');
    end
    
else %Does not exist
    disp('vgrid.in does not excist')
    nrvt=input('enter # of vertical layers ');
end

layerlist=cell(nrvt+1,1);
layerlist{1}='layer';
layerlist{2}='1 (Bottom)';
layerlist{end}=[num2str(nrvt) ' (Surface)'];
layerlist(3:end-1)= num2cell((2:nrvt-1));
layerlist61{1}='layer';
layerlist61{2}='Surface';

hpopupLayer = uicontrol('Style','popupmenu',...
    'String',layerlist,...
    'Position',[xDrop+sum(hxdrop(1:3)),yDrop,hxdrop(4),hy],...
    'Callback',{@popup_menu_Callback_layer});


hpopupAdd = uicontrol('Style','popupmenu',...
    'String',{'add','mesh','boundaries','coastline'},...
    'Position',[xDrop+sum(hxdrop(1:3))+hxdrop(4),yDrop,hxdrop(4),hy],...
    'Callback',{@popup_menu_Callback_add});

%add Zlevel window

Zcheck = uicontrol('Style','checkbox',...
    'String',  'Use Z level m :',...
    'Position',[xDrop,hxdrop(3)*10,hy*3/4,hxdrop(3)],...
    'Callback',{@Zcheckbox_Callback});

h = uicontrol('Style','edit','String','\default','Position',[xDrop,yDrop-hy*3/4+hxdrop(3)*5,hy*3/4,hxdrop(3)]);



ha = axes('Units','pixels','Position',[40,40,LX*0.9,LY*0.8]);
%align([hsurf,hmesh,hcontour,htext,hpopup],'Center','None');


% Initialize the GUI.
% Change units to normalized so components resize automatically.
%set([f,hplot,hpopup],'Units','normalized');
set(f,'Units','normalized');

ih=get(f,'children');

for i=1:length(ih)
    try
        ih(i)
        set(ih(i),'Units','normalized')
    catch msg
    end
end

%set(get(f,'children'),'Units','normalized')


% Generate the data to plot.

%% Global Variables
% Create a plot in the axes.
patchData = gr.hgrid.depth;


varname='depth';
prevvarname=' '; %previously selected varname 2 spaces
header=[];header2=[];
stack='stack';
time=0;
layer=1;
Dim=0;
prevDim=-1;

localData={}; %data sroage for scalar timeseries and slabs
localData2={}; %data sroage for v component timeseries and slabs if vector data
timeSeriesData={};
% Coordinates for time series
x=[];
y=[];
z=[];
%ph=[];%plot handle for location marker
gh=[];%handle to coast plotting
hq=[];%handle for quiver plotting
% Set positions for quiver makers

workSpaceExport=struct();

ph=gr_plot2(ll.hgrid,gr.hgrid.depth);
lh=[]; %legend handle
axe=axis; %spatial axis
akZoom();
% hold on
% if useShapefile
%     landplt=geoshow(shapefile);
%     set(findall(landplt,'type','patch'),'FaceColor',[.7 .7 .7]);
%     axis(axe)
% end

ch=colorbar;set(get(ch,'Ylabel'),'string',varname);
title(varname)

caxe=[]; %colorbar axis
[locs, dl]=arr_locs(axe,gr);
u=[];
v=[];
normal=[]; % absolute value sqrt(u^2+v^2)
u2=[];
v2=[];

%open and closed boundaries
openbd=[];
landbd=[];
toggle.plotbd=0;

%coastline Plot
Coast=[];
coastloaded=0;
toggle.plotcoast=0; % on/off switch for coast Plotting



% Menu variable for fortran extractions scripts
extractType=0;
%% callbacks
%  Pop-up menu callback. Read the pop-up menu Value property to
%  determine which item is currently displayed and make it the
%  current data. This callback automatically has access to
%  patchData because this function is nested at a lower level.



    function SetupButton_Callback(source,eventdata)
        
        fn.runDir = uigetdir([],'Select Selfe run directory - where hgrid.gr3 ');
        
        % open combined file
        fn.outDir = uigetdir(fn.runDir,'Select directory of combined binaries - eg 1_elev.61 ');
        
        gridName='hgrid.gr3';
        disp('loading grid file');gr.hgrid=gr_readHGrid(fullfile(fn.runDir,gridName));
        gridName='hgrid.ll';
        disp('loading grid file');ll.hgrid=gr_readHGrid(fullfile(fn.runDir,gridName));
        
        
        
        %varlist Different variables in folder
        varlist=dir([fn.outDir '/' stackList{2} '_*.6*' ]);
        varNames=cell(length(varlist)+1,1);
        varDims=zeros(length(varlist)+1,1);
        
        varNames{1}='depth';
        varDims(1)=0;
        
        for j=1:length(varlist);
            iunScore=find(varlist(j).name=='_');
            idot=find(varlist(j).name=='.');
            varNames{j+1}= varlist(j).name(iunScore+1:idot-1);
            varDims(j+1)=str2double(varlist(j).name(idot+1:end));
        end
        
        % count stacks
        %elevFiles=dir([fn.outDir '/*elev.61' ]);
        stackList=cell(length(dir([fn.outDir '/*elev.61' ]))+1,1);
        stackList{1}='stack';
        for i=1:length(stackList)-1;
            stackList{i+1}=num2str(i);
        end
        timesteplist={};
        timesteplist{1}='timestep';
        
        
        set(hpopup,'String',varNames);
        set(hpopupStack,'String',stackList)
        
        inregnodes=1:gr.hgrid.np;
        
        cla
        ph=gr_plot2(ll.hgrid,gr.hgrid.depth);
        axis([min(ll.hgrid.x) max(ll.hgrid.x) min(ll.hgrid.y) max(ll.hgrid.y)])
        axe=axis; %spatial axis
        akZoom();
        
        
    end


    function popup_menu_Callback(source,eventdata)
        
        % Determine the selected data set.
        str = get(source, 'String');
        val = get(source,'Value');
        % Set current data to the selected data set.
        
        
        switch str{val};
            case 'depth' % User selects Peaks.
                patchData = gr.hgrid.depth;
                layer=1;
                varname='depth';
                
                alternate=~alternate;
                if difference_plot && alternate
                    patchData = gr2.hgrid.depth-gr.hgrid.depth;
                end
                set(hpopupLayer,'String','')
                Plotbutton_Callback
            otherwise
                varname=varNames{val};
                Dim=varDims(val);
                
                
                %adjust layer selection menu
                if Dim > 62
                    set(hpopupLayer,'String',layerlist)
                    set(hpopupLayer,'Value',nrvt+1)
                    
                else
                    set(hpopupLayer,'String',layerlist61)
                    set(hpopupLayer,'Value',2)
                end
        end
        
        %reset stack, timestep and layer drop down
        set(hpopupStack,'value',1)
        set(hpopupTime ,'value',1)
        set(hpopupLayer ,'value',1)
    end




    function popup_menu_Callback_Stack(source,eventdata)
        % Determine the selected data set.
        %str = get(source, 'String');
        %val = get(source,'Value');
        str=get(hpopupStack,'String');
        val=get(hpopupStack,'Value');
        
        stack=str{val};
        
        %       end
        
        if ~strcmp(varname ,'depth') && ~strcmp(stack ,'stack')
            
            fileName=[stack '_' varname '.' num2str(Dim)];
            disp(['loading ' varname ' stack: ' stack ' (wait for printing done)'])
            
            %if ~strcmp(prevvarname,varname) % if variable changes or called first time fully load header
            % else only adapt path for  performane (data structure doesent change)
            if Dim~=prevDim;        % Do the same if Dimensions are euqal
                header=sz_readHeader(fullfile(fn.outDir,fileName));
            else
                header.fname=[header.fname(1:find(header.fname=='/',1,'last')) fileName];
            end
            %prevvarname=varname;
            prevDim=Dim;
            disp(['done loading ' varname ' stack: ' stack])
            
            %Load Timesteps
            if numel(timesteplist)==1
                for i=1:header.nSteps
                    timesteplist{i+1}=i;
                end
                set(hpopupTime,'String',timesteplist,'Value',1) % make timesteps available in dropdown
            end
        end
    end





    function popup_menu_Callback_Time(source,eventdata)
        if strcmp(stack,'stack')
            disp('please choose stack first')
            return
        end
        
        % Determine the selected data set.
        %str = get(source, 'String');
        %val = get(source,'Value');
        str = get(hpopupTime, 'String');
        val = get(hpopupTime,'Value');
        
        % Set current data to the selected data set.
        
        switch str{val};
            case 'timestep' % User selects Peaks.
                time=0;
                %
            otherwise
                
                time=str2double(str{val});
                
                if ~strcmp(varname ,'depth') && ~ strcmp(stack,'stack')
                    disp(['loading ' varname ' stack: ' stack ' timestep: ' num2str(time)])
                    [patchData, ts]=sz_readTimeStep(header,time); %
                    
                    patchData(patchData==dummy1)=nan;
                    patchData(patchData==dummy2)=nan;
                    
                    
                    if Dim==61
                        layer=1;
                        set(hpopupLayer,'value',2)
                    elseif Dim==62
                        layer=1;
                        u = patchData(:,1);%
                        v = patchData(:,2);%
                        patchData=sqrt(u.^2+v.^2);
                        set(hpopupLayer,'value',2)
                    elseif Dim==63
                        layer=nrvt;
                        [patchData] = map_sz2hts(header,patchData,1);%
                        set(hpopupLayer,'value',nrvt+1)
                    elseif Dim==64
                        layer=nrvt;
                        u = map_sz2hts(header,patchData(:,1),1);%
                        v = map_sz2hts(header,patchData(:,2),1);%
                        patchData=sqrt(u.^2+v.^2);
                        set(hpopupLayer,'value',nrvt+1)
                    end
                    
                    
                    if difference_plot
                        
                        fileName=header.fname(find(header.fname=='/',1,'last')+1:end);
                        header2=header;
                        header2.fname=fullfile(fn2.outDir,fileName);
                        
                        [temp ts]=sz_readTimeStep(header2,time); %
                        temp(temp==dummy1)=nan;
                        temp(temp==dummy2)=nan;
                        
                        if Dim==61
                            layer=1;
                            set(hpopupLayer,'value',2)
                        elseif Dim==62
                            layer=1;
                            u2 = temp(:,1);%
                            v2 = temp(:,2);%
                            u=u2-u;
                            v=v2-v;
                            temp=sqrt(u2.^2+v2.^2);
                            set(hpopupLayer,'value',2)
                        elseif Dim==63
                            layer=nrvt;
                            [temp] = map_sz2hts(header,temp,1);%
                            set(hpopupLayer,'value',nrvt+1)
                        elseif Dim==64
                            layer=nrvt;
                            u2 = map_sz2hts(header,temp(:,1),1);%
                            v2 = map_sz2hts(header,temp(:,2),1);%
                            u=u2-u;
                            v=v2-v;
                            temp=sqrt(u2.^2+v2.^2);
                            set(hpopupLayer,'value',nrvt+1)
                        end
                        
                        
                        patchData=temp-patchData;
                        
                    end
                    
                    
                    
                    
                    
                    Plotbutton_Callback
                    disp('done')
                end
                
                
                
        end
        
    end





    function popup_menu_Callback_add(source,eventdata)
        
        if strcmp(stack,'add')
            disp('please select what to add')
            return
        end
        
        % Determine the selected data set.
        str = get(source, 'String');
        val = get(source,'Value');
        % Set current data to the selected data set.
        
        switch str{val};
            case 'mesh' % User selects Peaks.
                axe=axis;
                addmesh(ll);
                axis(axe);
                
                if handles.tooltip.mesh==0
                    helpdlg('Klick inside an element to show its element and node numbers.')
                    handles.tooltip.mesh=1;
                end
                
                set (gcf, 'WindowButtonDownFcn', @mouseclick);
                
            case 'boundaries'
                
                toggle.plotbd=~toggle.plotbd;
                
                disp(['Add boundaries to plot = ' num2str(toggle.plotbd)])
                
                %plotBound2(gr.hgrid)
                if isempty(openbd)
                    [openbd, landbd]=plotBound(gr.hgrid);
                    PlotBD(openbd,landbd)
                else
                    PlotBD(openbd,landbd)
                end
            case 'coastline'
                
                toggle.plotcoast=~toggle.plotcoast; % toogle on off
                disp(['Add Landmask to plot = ' num2str(toggle.plotcoast)])
                
                if coastloaded==0
                    Coast=genCoast(gshhs_file);
                    coastloaded=1;
                    PlotCoast(Coast)
                else
                    if toggle.plotcoast
                        
                        
                        if ~isobject(gh)
                            PlotCoast(Coast)
                        else
                            set(gh,'visible','on')
                        end
                        
                    else
                        set(gh,'visible','off')
                    end
                    
                end
        end
        
        set(source,'Value',1);
        %set(source,'String',str{1});
        
    end



    function popup_menu_Callback_layer(source,eventdata)
        
        
        % Determine the selected data set.
        str = get(source, 'String');
        val = get(source,'Value');
        % Set current data to the selected data set.
        
        layer=val-1;
        
        if  time > 0 && ~strcmp(stack,'stack')
            Plotbutton_Callback
        else
            disp('First select a valid stack and time tep')
        end
        
        
    end



    function checkbox_Callback(source,eventdata)
        val = get(source,'Value');
        difference_plot=val;
        
        
        if val
            
            [gridName, fn2.runDir]=uigetfile([fn.runDir '/*.*'],'Select 2. hgrid file');
            
            
            % open combined file
            fn2.outDir = uigetdir(fn2.runDir,'Select 2. directory of combined binaries - eg 1_elev.61 ');
            
            
            gr2.hgrid=gr_readHGrid(fullfile(fn2.runDir,gridName));
            %gr2.hgrid=gr_readHGrid(fullfile(fn2.runDir,gridName));
            
            if(gr.hgrid.np~=gr2.hgrid.np || gr.hgrid.ne~=gr2.hgrid.ne)
                disp('grid dimensions of Setups to compare differ')
                set(source,'Value',0);
            end
        end
        video=[];
    end




    function VideoButton_Callback(source,eventdata)
        
        if axe~= axis
            axe=axis;
        end
        caxe=caxis;
        k=0;
        stack1=input('stack1 for video: ');
        stack2=input('stack2 for video: ');
        stack2=min(stack2,length(stackList)-1);
        
        if  isempty(video) || video.stack1~=stack1 || video.stack2~=stack2 || ...
                ~strcmp(video.varname,varname);
            video.varname=varname;
            video.stack1=stack1;
            video.stack2=stack2;
            for istack=stack1+1:stack2+1
                set(hpopupStack,'Value',istack)
                popup_menu_Callback_Stack
                for itime=2:length(timesteplist)
                    set(hpopupTime,'Value',itime)
                    popup_menu_Callback_Time
                    Plotbutton_Callback
                    caxis(caxe)
                    k=k+1;
                    pause(0.01)
                    video.patchData{k}=patchData(:,layer);
                    if Dim==62 || Dim == 64
                        video.u{k}=u(:,layer);
                        video.v{k}=v(:,layer);
                    end
                end
            end
            
        else % play again stored video
            if difference_plot
                
                rmsd=zeros(gr.hgrid.np,1);
            end
            
            for i=1:length(video.patchData)
                if difference_plot
                    rmsd=rmsd+video.patchData{i}.^2;
                end
                gr_plot2(ll.hgrid,video.patchData{i})
                caxis(caxe)
                
                
                if ~isempty(openbd)
                    %  hold on
                    for i=1:length(openbd)
                        plot(ll.hgrid.x(openbd{i}),ll.hgrid.y(openbd{i}));
                    end
                    for i=1:length(landbd)
                        plot(ll.hgrid.x(landbd{i}),ll.hgrid.y(landbd{i}),'r');
                    end
                    %    hold off
                end
                
                % Plot coast and Boundaries ?
                %if toggle.plotcoast;PlotCoast(Coast);end
                if toggle.plotbd;PlotBD(openbd,landbd);end;
                
                %                     if toggle.plotcoast
                %                         set(gh,'visible','on')
                %
                %                     else
                %                         set(gh,'visible','off')
                %                     end
                
                
                if Dim==62 || Dim == 64
                    %                     if isobject(qh)
                    %                     delete(qh)
                    %                     end
                    %
                    if isgraphics(hq);delete(hq);end
                    %hold on
                    qh=quiver(ll.hgrid.x(locs),ll.hgrid.y(locs),video.u{i}(locs),video.v{i}(locs),'k');
                    %hold off
                end
                
                
                pause(0.01)
            end
            if difference_plot
                rmsd=sqrt(rmsd/length(video.patchData));
                figure(1)
                gr_plot2(ll.hgrid,rmsd)
                
                % Plot coast and Boundaries ?
                if toggle.plotcoast;PlotCoast(Coast);end
                if toggle.plotbd;PlotBD(openbd,landbd);end;
                
                title('rmsd')
                %video.rmsd=rmsd;
                assignin('base',[varname '_rmsd'],rmsd)
            end
        end
        
    end

    function Plotbutton_Callback(source,eventdata)
        
        if axe~=axis
            axe=axis;
            
            %set limits to zoomed area
            inregnodes=ll.hgrid.x >= axe(1) & ll.hgrid.x <= axe(2) & ...
                ll.hgrid.y >= axe(3) & ll.hgrid.y <= axe(4);
            %else
            %   inregnodes(:)=true;
        end
        caxisLim=[min(patchData(inregnodes,layer)) max(patchData(inregnodes,layer))];
        
        
        % Display surf plot of the currently selected data.
        if ~strcmp(varname,'depth')
            if strcmp(stack,'stack')
                disp('please firt select stack')
                return
            end
            if ~time
                disp('please firt select timestep')
                return
            end
        end
        
        %cla
        workSpaceExport.patch=patchData(:,layer);
        assignin('base','result',workSpaceExport)
        %        if isgraphics(ph);delete(ph);end
        if isobject(ph);delete(ph);end
        if ~isempty(lh);delete(lh);end
        
        ph=gr_plot2(ll.hgrid,patchData(:,layer));
        ch=colorbar;set(get(ch,'Ylabel'),'string',varname);
        axis(axe)
        
        
        
        % Plot coast and Boundaries ?
        if toggle.plotcoast;PlotCoast(Coast);end
        if toggle.plotbd;PlotBD(openbd,landbd);end;
        
        if difference_plot
            if ~strcmp(varname,'depth') || alternate
                title(['\Delta' varname])
                colormap redblue
            else
                title([varname ' 1'])
                colormap jet
            end
            caxis([-1 1]*max(abs(patchData(:,layer))) )
        elseif ~strcmp(varname,'depth')
            datum=datenum(header.startTime,'mm/dd/yyyy hh:MM:ss')+((str2double(stack)-1)*header.nSteps+time)*header.dt/86400;
            title([varname ' @' datestr(datum,'mm/dd/yyyy hh:MM:ss')]);
        end
        
        if Dim==62 || Dim==64
            
            [locs, dl]=arr_locs(axe,gr);
            
            
            hold on
            if ~difference_plot
                normal=patchData(:,layer);
                normal(normal==0)=1;
                %if isgraphics(hq);delete(hq);end
                
                %if ~isempty(hq);delete(hq);end
                if isobject(hq);delete(hq);end
                hq=quiver(ll.hgrid.x(locs),ll.hgrid.y(locs),u(locs,layer)./normal(locs),v(locs,layer)./normal(locs),1,'k');
            else
                maxx=xlim;
                maxx=maxx(2)*0.9;
                maxy=ylim;
                maxy=(maxy(2)-maxy(1))*.1+maxy(1);
                
                if varname=='wind'
                    uval=10;
                else
                    uval=1;
                end
                
                %if isgraphics(hq);delete(hq);end
                if ~isempty(hq);delete(hq);end
                
                hq=quiver([ll.hgrid.x(locs); maxx],[ll.hgrid.y(locs); maxy],[u(locs,layer); uval],[v(locs,layer); 0],1,'k');
                text(maxx,maxy*1.002,sprintf('%i m/s',uval))
            end
            
          
            
        end
        
            %transect exists
            if isfield(transect,'nor')
                
            end
        
        caxis(caxisLim);
    end


    function ExtractButton_Callback(source,eventdata)
        % Display surf plot of the currently selected data.
        
        if   strcmp(varname,'depth')
            disp('please select variable other than depth')
            return
        end
        
        datadir=fn.outDir;
        binary_file=[varname '.' num2str(Dim)];
        [extractType, XYZ, ics, ifs, istack1, istack2, varargout]=invokeExtract(fn,gr,ll,binary_file,datadir,days);
        
         hovmue=varargout;
        
        if isempty(hovmue) || strcmp(hovmue.type,'transect')
            
            [d,d]=initialize_fortran_extract(fn,gr,ll,binary_file,ifs,istack1,istack2,XYZ);
            !mv SubRoutines/fortran/ts_extract/fort.18 SubRoutines/fortran/ts_extract/fort.18_1
            !mv SubRoutines/fortran/ts_extract/fort.19 SubRoutines/fortran/ts_extract/fort.19_1
            !mv SubRoutines/fortran/ts_extract/fort.20 SubRoutines/fortran/ts_extract/fort.20_1
            
            if difference_plot
                
                [d,d]=initialize_fortran_extract(fn2,gr,ll,binary_file,ifs,istack1,istack2,XYZ);
                !mv SubRoutines/fortran/ts_extract/fort.18 SubRoutines/fortran/ts_extract/fort.18_2
                !mv SubRoutines/fortran/ts_extract/fort.19 SubRoutines/fortran/ts_extract/fort.19_2
                !mv SubRoutines/fortran/ts_extract/fort.20 SubRoutines/fortran/ts_extract/fort.20_2
                
                
                
            end
            
        end
        
        if  isempty(varargout) % no hovmueller
            plot_extract
        else
            
            if strcmp(hovmue.type,'transect')
                plot_hovmueller(ll,hovmue,dummy1,dummy2,Dim,nrvt,varname)
            else
                profile on
                hovmuellerProfile(fn,ll,hovmue,dummy1,dummy2,Dim,nrvt,varname);
                profile viewer
                profile off
            end
            
        end
        
    end




    function CreateButton_Callback(source,eventdata)
        % Display surf plot of the currently selected data.
        
        
        
        
        extractType=menu('extract options','transect');
        
        if extractType==1
            
            if strcmp(varname,'depth')
                disp('please choose a variable other than depth')
                
            elseif ~(Dim==64 || Dim==63)
                disp('please choose a variable with vertical extend')
            else
                transect=TranSectProfile(fn,ll,gr,[varname '.' num2str(Dim)]);
                assignin('base','transect',transect)
            end
        end
        
        
    end


    function exportButton_Callback(source,eventdata)
        caxe=caxis;
        axe=axis;
        options={'Plot'};
        if ~isempty(video)
            options{end+1}='video';
        end
        if ~isempty(timeSeriesData)
            options{end+1}='timeseries';
        end
        
        
        select=menu('export',options);
        if select==1
            
            
            gcf=figure(1);
            h = findobj(gcf,'type','axes'); % Find the axes object in the GUI
            f1 = figure; % Open a new figure with handle f1
            s = copyobj(h(2),f1); % Copy axes object h into figure f1
            set(s,'Units','normalized','Position',[.1 .05 .85 .9])
            ch=colorbar;set(get(ch,'Ylabel'),'string',varname);
            %
            %             figure
            %             gr_plot(gr.hgrid,patchData)
            %             hold on
            %               if Dim==62 || Dim==64
            %                   hold on
            %                   %normal=patchData(:,layer);
            %                   quiver(gr.hgrid.x(locs),gr.hgrid.y(locs),u(locs)./patchData(locs,layer),...
            %                       v(locs)./patchData(locs,layer),1,'color','k')
            %                   colorbar
            %               end
            %
            %               if ~isempty(openbd)
            %                 keyboard
            %               end
            
        else %video_export
            axe=axis;
            video.gr=gr;
            disp('exporting video')
            videoDir=uigetdir('select_folder');
            %writerObj =VideoWriter(fullfile(videoDir,[varname '_' num2str(video.stack1) '_' num2str(video.stack1) '.avi']));
            %set(writerObj,'FrameRate',15);
            %set(writerObj,'Quality',100);
            
            dt= header.dt;
            nSteps=header.nSteps;
            startdate=datenum(header.startTime,'mm/dd/YYYY HH:MM:SS');
            startdate=startdate+(video.stack1-1)*nSteps*dt/86400;
            
            n=length(video.patchData);
            fh=figure;
            input('set to full screen and press return ')
            set(fh,'Visible','off');
            %open(writerObj);
            increment=dt/86400;
            datum=startdate;
            
            do_quiv=0;
            if Dim==62 || Dim==64
                %if axe~=axis
                %   axe=axis;
                [locs, dl]=arr_locs(axe,gr);
                do_quiv=1;
                %end
            end
            
            if difference_plot
                   colormap redblue 
            end
            
            for i=1:n
                %if isgraphics(ph)
                if isobject(ph)
                    delete(ph)
                end
                if isobject(lh)
                    delete(lh)
                end
                
                ph=gr_plot2(ll.hgrid,video.patchData{i});
                hold on
                
                % Plot coast and Boundaries ?
                %if toggle.plotcoast;PlotCoast(Coast);end
                if toggle.plotbd;PlotBD(openbd,landbd);end;
                
                %Boundary
                if ~isempty(openbd)
                    PlotBD(openbd,landbd)
                end
                
                title(datestr(datum))
                caxis(caxe)
                ch=colorbar;set(get(ch,'Ylabel'),'string',varname);
                fname=sprintf('%03i_%s',i,varname);
                axis(axe)
                
                
                
                
                if do_quiv;
                    %if isgraphics(hq)
                    if isobject(hq)
                        delete(hq)
                    end
                    hold on
                    hq=quiver(ll.hgrid.x(locs),ll.hgrid.y(locs),video.u{i}(locs),video.v{i}(locs),'k');
                end
                
                
                
                
                
                print(fh,'-dpng',fullfile(videoDir,fname));
                %if isgraphics(ph);delete(ph);end
                if isobject(ph);delete(ph);end
                if do_quiv;
                    %if isgraphics(hq);delete(hq);end
                    if isobject(hq);delete(hq);end
                end
                %    frame = getframe(fh);
                %     writeVideo(writerObj,frame);
                datum=datum+increment;
            end
            close(fh);
            % close(writerObj);
            save(fullfile(videoDir,'video.mat' ),'video');
            disp('done xporting video')
            
        end
        
        
        
    end





%% Other functions





%quiver arrow origins
    function [locs, dl]=arr_locs(axe,gr)
        
        dl=(axe(2)-axe(1))/50;
        
        locs=zeros(length(axe(1):dl:axe(2))*length(axe(3):dl:axe(4)),1);
        k=1;
        for lon=axe(1):dl:axe(2)
            
            for lat=axe(3):dl:axe(4)
                
                %nearest grid cell
                [was, wo]=min(abs(ll.hgrid.x-lon) + abs(ll.hgrid.y-lat) );
                locs(k)=wo;
                k=k+1;
            end
            
            %  icrit=find(gr.hgrid.x >=axe(1) & gr.hgrid.x <=axe(2) & gr.hgrid.y >=axe(3) & gr.hgrid.y <=axe(4) ...
            %     & c(:,end) > 4);
        end
        
    end



%generate Landmask from binary gshhs files using matlab mapping toolbox
    function Coast=genCoast(gshhs_file)
        
        disp('Extracting coastline polygons from gshhs database - this might take a while')
        latlim=axe(3:4);
        lonlim=axe(1:2);
        Coast = gshhs(gshhs_file, latlim, lonlim);
        
        %reduce Polygon extend for faster Plotting
        for i=1:length(Coast)
            
            [Coast(i).Lat,Coast(i).Lon] = maptrimp(Coast(i).Lat,Coast(i).Lon,latlim + [-1 1],lonlim + [-1 1]);
        end
    end

    function PlotCoast(Coast)
        
        %if isgraphics(gh)
        if isobject(gh)
            
        else
            axis(axe)
            hold on
            gh=geoshow(Coast);
            hold off
            axis(axe)
            
            
        end
    end

    function PlotBD(openbd,landbd)
        hold on
        for j=1:length(openbd)
            hb1(j)=plot(ll.hgrid.x(openbd{j}),ll.hgrid.y(openbd{j}),'linewidth',1);
        end
        for j=1:length(landbd)
            hb2(j)=plot(ll.hgrid.x(landbd{j}),ll.hgrid.y(landbd{j}),'k','linewidth',1);
            
        end
        
        hg1 = hggroup;
        set(hb1,'Parent',hg1)
        set(hg1,'Displayname','Open Boundaries')
        hg2 = hggroup;
        set(hb2,'Parent',hg2)
        set(hg2,'Displayname','Closed Boundaries')
        lh=legend([hg1 hg2]);
        hold off
        
        
        
    end




    function plot_extract
        
        switch extractType
            
            case 1%'extract time series
                localData={};
                
                localData{1}=load('fort.18_1');
                localData{1}(localData{1}==dummy1)=nan;
                localData{1}(localData{1}==dummy2)=nan;
                %tstime=localData(:,1);
                % localData{1}=localData{1}(:,[1 2:nrvt:end]);
                
                % assign the value of the dot product to the variable 'result' in the workspace
                timeSeriesData{1}=localData{1}(:,2:end);
                
                if difference_plot
                    localData{2}=load('fort.18_2');
                    localData{2}(localData{2}==dummy1)=nan;
                    localData{2}(localData{2}==dummy2)=nan;
                    %tstime=localData(:,1);
                    %localData{2}=localData{2}(:,[1 2:nrvt:end]);
                    % assign the value of the dot product to the variable 'result' in the workspace
                    timeSeriesData{2}=localData{2}(:,2:end);
                end
                
                
                figure(1)
                if Dim==62 || Dim==64
                    
                    if (Dim==62)
                        nrvti=1;
                    else
                        nrvti=nrvt;
                    end
                    
                    
                    localData2{1}=load('fort.19_1');
                    localData2{1}(localData2{1}==dummy1)=nan;
                    localData2{1}(localData2{1}==dummy2)=nan;
                    localData2{1}=localData2{1}(:,[1 2:nrvti:end]);
                    timeSeriesData{1}=timeSeriesData{1}+1i*localData2{1}(:,2:end);
                    
                    if difference_plot
                        localData2{2}=load('fort.19_2');
                        localData2{2}(localData2{2}==dummy1)=nan;
                        localData2{2}(localData2{2}==dummy2)=nan;
                        localData2{2}=localData2{2}(:,[1 2:nrvti:end]);
                        timeSeriesData{2}=timeSeriesData{2}+1i*localData2{2}(:,2:end);
                    end
                    
                    
                    figure(3)
                    clf
                    LegendCell={};
                    for i=1:1+difference_plot
                        
                        if difference_plot
                            subplot(1,2,i)
                        end
                        plot(timeSeriesData{i})
                        hold on
                        plot(timeSeriesData{i}(1,:),'o');
                        plot(timeSeriesData{i}(end,:),'x');
                        xlabel('u')
                        ylabel('v')
                        
                        for ii=size(LegendCell,1)+1:size(timeSeriesData,2)*i
                            LegendCell{ii}=sprintf('P %i',ii);
                        end
                        
                        legend(LegendCell)
                    end
                    
                end
                
                %make Data available in Workspace
                varname(varname=='-')='_';
                eval(['workSpaceExport.' varname '=timeSeriesData' ]);
                %assignin('base','result',{localData,localData})
                assignin('base','result',workSpaceExport)
                
                %  if i==(1+difference_plot)
                figure(2)
                clf
                hold on
                %plot(localData(:,1))
                
                if ~difference_plot
                    plot(localData{1}(:,2:end),'.-')
                    
                    LegendCell={};
                    ientry=0;
                    
                    if Dim==62 || Dim==64
                        plot(localData2{1}(:,2:end),'--')
                        
                        for ii=1:size(timeSeriesData{1},2)
                            ientry=ientry+1;
                            LegendCell{ientry}=sprintf('P %i u',ii);
                        end
                        
                        for ii=1:size(timeSeriesData{1},2)
                            ientry=ientry+1;
                            LegendCell{ientry}=sprintf('P %i v',ii);
                        end
                        
                        %legend('u','v');
                        %linkaxes([ax1,ax2],'x')
                    else
                        
                        ientry=0;
                        for ii=1:size(timeSeriesData{1},2)
                            ientry=ientry+1;           end
                        %legend(varname)
                        
                    end
                    
                    %linkaxes([ax1,ax2],'x')
                    legend(LegendCell)
                    
                else
                    
                    if Dim==62 || Dim==64
                        
                        ax1=subplot(2,2,1);
                        plot(localData{1}(:,2:end),'.-')
                        hold on
                        plot(localData{2}(:,2:end),'--o')
                        title([varname 'u'])
                        
                        LegendCell={};
                        ientry=0;
                        addendum=['','_setup2'];
                        for i=1:size(timeSeriesData,2)
                            for ii=1:size(timeSeriesData{1},2)
                                ientry=ientry+1;
                                LegendCell{ientry}=sprintf('u @ P%s %i',addendum(i),ientry);
                            end
                        end
                        legend(LegendCell)
                        
                        ax2=subplot(2,2,2);
                        plot(localData2{1}(:,2:end),'.-')
                        hold on
                        plot(localData2{2}(:,2:end),'--o')
                        title([varname 'v'])
                        
                        LegendCell={};
                        ientry=0;
                        addendum=['','_setup2'];
                        for i=1:size(timeSeriesData,2)
                            for ii=1:size(timeSeriesData{1},2)
                                ientry=ientry+1;
                                LegendCell{ientry}=sprintf('v @ P%s %i',addendum(i),ientry);
                            end
                        end
                        legend(LegendCell)
                        
                        
                        
                        LegendCell={};
                        ax3=subplot(2,2,3);
                        plot(localData{1}(:,2:end)-localData{2}(:,2:end))
                        for ii=1:size(timeSeriesData{1},2)
                            LegendCell{ii}=sprintf('Delta u @P %i',ii);
                        end
                        legend(LegendCell)
                        title(['\Delta ' varname])
                        %legend(varname)
                        
                        
                        LegendCell={};
                        ax4=subplot(2,2,4);
                        plot(localData2{1}(:,2:end)-localData2{2}(:,2:end))
                        for ii=1:size(timeSeriesData{1},2)
                            LegendCell{ii}=sprintf('Delta v @P %i',ii);
                        end
                        legend(LegendCell)
                        title(['\Delta ' varname])
                        %legend(varname)
                        
                        linkaxes([ax1,ax3],'x')
                        linkaxes([ax2,ax4],'x')
                        
                    else
                        
                        ax1=subplot(2,1,1);
                        plot(localData{1}(:,2:end),'.-')
                        hold on
                        plot(localData{2}(:,2:end),'--o')
                        title(varname)
                        
                        LegendCell={};
                        ientry=0;
                        addendum=['','_setup2'];
                        for i=1:size(timeSeriesData,2)
                            for ii=1:size(timeSeriesData{1},2)
                                ientry=ientry+1;
                                LegendCell{ientry}=sprintf('P%s %i',addendum(i),ientry);
                            end
                        end
                        legend(LegendCell)
                        LegendCell={};
                        ax2=subplot(2,1,2);keyboard
                        plot(localData{1}(:,2:end)-localData{2}(:,2:end))
                        for ii=1:size(timeSeriesData{1},2)
                            LegendCell{ii}=sprintf('Delta P %i',ii);
                        end
                        legend(LegendCell)
                        title(['\Delta ' varname])
                        %legend(varname)
                        linkaxes([ax1,ax2],'x')
                    end
                    
                end
                
                
                
            case 2%'extract S-level Layer
                
                
                %vectordata
                if Dim == 62 || Dim == 62
                    
                    localData=load(...
                        sprintf('Subroutines/fortran/%i_%s.slab%i',...
                        9,varname,41));
                    
                    localData(localData==dummy2 | localData==dummy1)=nan;
                    figure
                    gr_plot2(ll.hgrid,localData(:,end))
                    
                    
                else %Scalar
                    %Dimensions [gridpoints timesteps]
                    
                    
                end
                
                
            case 3%'extract Z level Layer
                
                keyboard
                
            case 4% 'extract values at xyt
                
                keyboard
                
            case 5%extract values at xyzt'
                
                keybord
        end
        
        
    end



%keypress Function for figure
    function [NowSelected] = keypress(src,evnt)
        
        %eneable keboard controol to scroll through timesteps
        %by using key arrows
        
        if ~strcmp(varname,'depth')
            
            
            if strcmp(evnt.Key,'rightarrow')
                % was du willst #1
                if get(hpopupStack,'Value')==1  %%initial first timestep/stack
                    set (hpopupStack,'Value',2)
                    popup_menu_Callback_Stack(src,evnt)
                    set(hpopupTime,'Value',2)
                    popup_menu_Callback_Time(src,evnt)
                elseif time<header.nSteps  %before last timestep
                    set(hpopupTime,'Value',get(hpopupTime,'Value')+1)
                    popup_menu_Callback_Time(src,evnt)
                    
                elseif str2double(stack)<str2double(stackList(end)) %before last stack
                    set (hpopupStack,'Value',get(hpopupStack,'Value')+1)
                    popup_menu_Callback_Stack(src,evnt)
                    set(hpopupTime,'Value',2)
                    popup_menu_Callback_Time(src,evnt)
                else %begin fromstart
                    set (hpopupStack,'Value',2)
                    popup_menu_Callback_Stack(src,evnt)
                    set(hpopupTime,'Value',2)
                    popup_menu_Callback_Time(src,evnt)
                end
                
            end
            % was du willst #1
            
            
            
            %going backward
        elseif strcmp(evnt.Key,'lefttarrow')
            
            
            if get(hpopupStack,'Value')==1   %initial first timestep/stack
                set (hpopupStack,'Value',length(stackList))
                popup_menu_Callback_Stack(src,evnt)
                set(hpopupTime,'Value',length(timesteplist))
                
                popup_menu_Callback_Time(src,evnt)
            elseif time>1  %after first timestep
                set(hpopupTime,'Value',get(hpopupTime,'Value')-1)
                popup_menu_Callback_Time(src,evnt)
                
            elseif str2double(stack)>str2double(stackList(2)) %greater first stack
                set (hpopupStack,'Value',get(hpopupStack,'Value')-1)
                popup_menu_Callback_Stack(src,evnt)
                set(hpopupTime,'Value',length(stackList))
                popup_menu_Callback_Time(src,evnt)
            else %begin from end
                set (hpopupStack,'Value',length(stackList))
                popup_menu_Callback_Stack(src,evnt)
                set(hpopupTime,'Value',length(timesteplist))
            end
            
            
        end
        
    end


% Mouseover
%set (gcf, 'WindowButtonMotionFcn', @mouseover);

    function [data] = mouseclick(gcbo,eventdata,handles)
        %function [data] = mouseover()
        
        c = get (gca, 'CurrentPoint'); % get mouse coordinates
        %[val pos]=min(abs(ll.hgrid.nodes(:,2)-c(2,1))+abs(ll.hgrid.nodes(:,3)-c(2,2)));
        
        parents=find_parent(ll,c(2,:));
        
        hold on
        if ll.hgrid.elem(parents,2)==3
            
            elmNodes=ll.hgrid.elem(parents,3:5);
            ex=mean(ll.hgrid.x(ll.hgrid.elem(parents,3:5)));
            ey=mean(ll.hgrid.y(ll.hgrid.elem(parents,3:5)));
            
        elseif ll.hgrid.elem(parents,2)==4
            
            elmNodes=ll.hgrid.elem(parents,3:6);
            ex=mean(ll.hgrid.x(ll.hgrid.elem(parents,3:6)));
            ey=mean(ll.hgrid.y(ll.hgrid.elem(parents,3:6)));
            
            
        end
        plot(ll.hgrid.x(elmNodes([1:end 1])),ll.hgrid.y(elmNodes([1:end 1])),'r--');
        plot(ex,ey,'r.');
        text(ex,ey,sprintf('element: %i',parents),'fontsize',18);
        plot(ll.hgrid.x(elmNodes),ll.hgrid.y(elmNodes),'ro');
        text(ll.hgrid.x(elmNodes),ll.hgrid.y(elmNodes),num2str((elmNodes)','Node: %i'),'fontsize',18);
        hold off
        
    end

end
