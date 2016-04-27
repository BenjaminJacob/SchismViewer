

function invokeExtract_Sslab(gr,ll,fill_in,basedir,binary_file,ialong_S,istack1,istack2,klev0)


%% get XYZ(points)



%% create extract optionsFile for fortran routine
fid=fopen('SubRoutines/fortran/extract_options','w');
fprintf(fid,'fill_in=%i       !Fill in value used for abnormal case (e.g., below bottom, dry etc)\n',fill_in)
fprintf(fid,'basedir="%s"   !basedir - directory where binary and build point files are stored\n',basedir)
fprintf(fid,'binary_file=%s !binary_file - elev.61, hvel.64 etc; it should be inside basedir\n',binary_file)


fprintf(fid,'ialong_S=%i ! ialng_S 1) 0:= along Zelvel 1) 1:= along Slevel  ;\n',ialong_S)
if ialong_S==0
    zout=input('Enter Zlevel to extract slab for : ');
    fprintf(fid,'zout=%f ! Zlevel of extracted slab \n',zout)
    fprintf(fid,'ifs=%i !uare input z coordinates 0) z cordinates or 1) relative to free surface; and in this case, make sure all z>=0 ? Enter 0 or 1: ',ifs);
    ifs=input('are input z coordinates 0) z cordinates or 1) relative to free surface; and in this case, make sure all z>=0 ? Enter 0 or 1: ');
    fprintf(fid,'ifs=%i !uare input z coordinates 0) z cordinates or 1) relative to free surface; and in this case, make sure all z>=0 ? Enter 0 or 1: ',ifs);
else
   % klev0=input('Enter S-level to extract slab for : ');
    fprintf(fid,'klev0=%i ! Slevel of extracted slab\n',klev0)
end


fprintf(fid,'istack1=%i !     istack[1,2] - start and end stack # for binary outputs;\n',istack1)
fprintf(fid,'istack2=%i\n',istack2)

end