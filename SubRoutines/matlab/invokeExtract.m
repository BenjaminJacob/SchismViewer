
%% SElect

function [extractType XYZ ics ifs istack1 istack2 varargout]=invokeExtract(fn,gr,ll,varnames,datadir,days,varargin)



extractType=menu('extract options','extract time series','extract S-level Layer','extract Z level Layer','extract transect profile','Hovmöller diagram');


%keyboard
if nargin>4 && extractType==2
    varnames={varnames 'zcor.63'};
else
    varnames={varnames};
end


for ivarname=1:length(varnames)
    
    varname=varnames{ivarname};
    fill_in=-9999; % Fill value for empty data
    
    
    basedir=datadir;%directory of combined ressults
    binary_file=varname;%variable name with extension e.g. elev.61
    
    
    switch extractType
        
        case 1%'extract time series
            
            %bp_file='station.ts';
            bp_file='station.bp';
            bp_fileName='bp_file';  % Build PointFIleName
            
            [XYZ ics ifs istack1 istack2]=invokeExtract_TS3(fn,gr,ll,fill_in,basedir,binary_file,bp_fileName,bp_file,days);
            prog='read_output7M_xyz';
            varargout{1}=[];
            
            
            
            
        case 2%'extract S-level Layer
            
            ialong_S=1;
            klev0=input('Enter S-level to extract slab for : ');
            istack1=input('enter stack 1: ');
            istack2=input('enter stack 2: ');
            invokeExtract_sslab(gr,ll,fill_in,basedir,binary_file,ialong_S,istack1,istack2,klev0)
            prog='extract_slab';
            
        case 3%'extract Z level Layer
            
            keyboard
            
            ialong_S=0;
            invokeExtract_slab(gr,ll,fill_in,basedir,binary_file,ialong_S)
            prog='extract_slab';
        case 4% 'extract values at xyt
            
            transect=TranSectProfile(fn,ll,gr,varname);
            assignin('base','transect',transect)
            
            
        case 5%extract values at xyzt'
            
            
            hovmue=ExtractHovmueller(fn,gr,ll,varname,days);
 
            if strcmp(hovmue.type,'transect')
                prog='read_output7M_xyz';
                bp_file='station.bp';
                bp_fileName='bp_file';  % Build PointFIleName
                XYZ=hovmue.XYZ;
                ics=hovmue.ics;
                ifs=hovmue.ifs;
                
                
            else
                ics=0;
                ifs=0;
            end
            XYZ=0;
            istack1=hovmue.istack1;
            istack2=hovmue.istack2;
            varargout{1}=hovmue;
    end
    
    %% RUN FORTRAN SCRIPTS
    
    %     %temporary change matlabs librry path to system path to get script work
    %     matlbPath=getenv('LD_LIBRARY_PATH');
    %     %matlbPath='/programme/matlab/sys/os/glnxa64:/programme/matlab/bin/glnxa64:/programme/matlab/extern/lib/glnxa64:/programme/matlab/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/programme/matlab/sys/java/jre/glnxa64/jre/lib/amd64/server:/programme/matlab/sys/java/jre/glnxa64/jre/lib/amd64'
    %
    %     %replace matlablibary with system library to get gcc compiled fortran
    %     %scripts work
    %     setenv('LD_LIBRARY_PATH','/usr/lib')
    %
    %     maindir=pwd;
    %     try
    
    %         %execute fortran Program
    %         cd('SubRoutines/fortran');
    %         system(sprintf('./%s',prog))
    %
    %     catch
    %         setenv('LD_LIBRARY_PATH',matlbPath)
    %         cd maindir
    %     end
    %
    %     cd(maindir)
    %     %reset to default
    %     setenv('LD_LIBRARY_PATH',matlbPath)
    
end % varname loop




