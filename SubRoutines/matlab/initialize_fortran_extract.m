function [timevec ts_extract]=initialize_fortran_extract(fn,gr,ll,binary_file,ifs,stack1,stack2,XYZ)

% function to initialize and call fortran extraction routine
% input:
%1) fn  - fn.runDir:= >directory where setup resides <
%       - fn.outDir:= >directory whre combined binaries reside (e.g. elev.61)<
%2) gr  - gridstructure created with gr.hgrid=gr_readHgrid from (hgrid.gr3)
%3) ll  - gridstructure created with ll.hgrid=gr_readHgrid from (hgrid.ll)
%4) binary_file - name of binaryfile to extract e.g (elev.61,salt.63, etc)
%5) stack1 - beginning stack for extraing timeseries
%6) stack2 - ending stack for extraing timeseries
%7) ifs - specification of Z-value type ifs=1: Zcoordinates  ifs=2: %coordinates relative to free surface
%  XYZ  - coordinate matrix of dimension [number of points, cordinate of point (xyz)]

%cartesian or lonlat
if gr.hgrid.x(1)~=ll.hgrid.x(1)
    ics=1;
else
    ics=2;
end

fill_in=-9999;
bp_file='station.bp';
npoints=size(XYZ,1);

%% extract options
fid=fopen(fullfile('SubRoutines/fortran/ts_extract','extract_options'),'w');
fprintf(fid,'fill_in=%i       !Fill in value used for abnormal case (e.g., below bottom, dry etc)\n',fill_in);
fprintf(fid,'rundir="%s"   !directory where grid files are stored (e.ge one level above outputs)\n',fn.runDir);
fprintf(fid,'basedir="%s"   !basedir - directory where binary and build point files are stored\n',fn.outDir);
fprintf(fid,'binary_file=%s !binary_file - elev.61, hvel.64 etc; it should be inside basedir\n',binary_file);
fprintf(fid,'istack1=%i !     istack[1,2] - start and end stack # for binary outputs;\n',stack1);
fprintf(fid,'istack2=%i\n',stack2);
fprintf(fid,'%s=%s !  %s - build point file (see station.ts.sample / station.xyzt.sample)\n','bp_file',bp_file,'bp_file');
fprintf(fid,'%s=%i  ! (ics=1: Cartesian coordinates; =2: lat/lon.\n','ics',ics);
fprintf(fid,'%s=%i !ifs=0: input z are z coordinates;ifs =1: input z are relative to free surface; and in this case, make sure all z>=0\n','ifs',ifs);
fclose(fid)
%% create station.ts input file for fortran routine


itransect=1; %0 time seires;  1: transect %funktionier nur mit 1
fid=fopen(fullfile('SubRoutines/fortran/ts_extract',bp_file),'w');
%fprintf(fid,'%i %i %i  !itransect,ics,ifs (itransect=0: time series; =1:transect. ics=1: Cartesian coordinates; =2: lat/lon. ifs=0: input z are z coordinates; =1: input z are relative to free surface; and in this case, make sure all z>=0 :\n',itransect,ics,ifs);
fprintf(fid,'%i !total # of points\n',npoints);
fprintf(fid,'%i %14.6f %14.6f %4.4f !point #, x, y, z\n', XYZ(1,:));
if npoints>1
    fprintf(fid,'%i %14.6f %14.6f %4.4f\n', XYZ(2:end,:)');
end
fclose(fid);

%run fortran code

%% RUN FORTRAN SCRIPTS

%temporary change matlabs librry path to system path to get script work
matlbPath=getenv('LD_LIBRARY_PATH');
%matlbPath='/programme/matlab/sys/os/glnxa64:/programme/matlab/bin/glnxa64:/programme/matlab/extern/lib/glnxa64:/programme/matlab/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/programme/matlab/sys/java/jre/glnxa64/jre/lib/amd64/server:/programme/matlab/sys/java/jre/glnxa64/jre/lib/amd64'

%replace matlablibary with system library to get gcc compiled fortran
%scripts work
setenv('LD_LIBRARY_PATH','/usr/lib')

maindir=pwd;
try

    %execute fortran Program
    cd('SubRoutines/fortran/ts_extract')
    !./read_output7M_xyz
    cd(maindir);
catch
    setenv('LD_LIBRARY_PATH',matlbPath)
    cd(maindir)
end

cd(maindir)
%reset to default
setenv('LD_LIBRARY_PATH',matlbPath)


ts_extract=load('fort.18');
timevec=ts_extract(:,1);
ts_extract=ts_extract(:,2:end);