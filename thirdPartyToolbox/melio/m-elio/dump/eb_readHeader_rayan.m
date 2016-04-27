function [ header ] = eb_readHeader_new( fname , modelflag)
% [ header ] = eb_readHeader( fname )
% header for the elcirc binary fname
% modelflag = 1: SELFE
% modelflag = 0: ELCIRC
%
% Sergey Frolov mar 08, 2004
% Modified Ryan W. Kilgren Mar 30, 2005

% tic
header.fname        = fname;
header.idx.flag     = 0;
header.stepSize     = 0;

fid                 = fopen(fname);
header.fid          = fid;
header              = readHHeader(fid,header,modelflag);
% if(modelflag == 0)
% header              = readVgrid_elcirc(fid,header,modelflag);
% elseif(modelflag == 1)
% header              = readVgrid_selfe(fid,header,modelflag);
% end;
header              = readHgrid(fid,header,modelflag);
header.dataStartPos = ftell(fid);
header.fid          = fclose(fid);

% toc

function [header]=readHHeader(fid,header,modelflag)
	header.dataFormat   = char(fread(fid,48,'uchar')');
	header.version      = char(fread(fid,48,'uchar')');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Trying to make this function more dynamic. %
    % If you fix this if statement I think it    %
    % should work dynamically.                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     if(length(findstr(header.version,'ELCIRC')))
%         modelflag = 0;
%     elseif(length(findstr(header.version,'SELFE')))
%         modelflag = 1;
%     end
	header.startTime    = char(fread(fid,48,'uchar')');
	header.varType      = char(fread(fid,48,'uchar')');
	header.varDimension = char(fread(fid,48,'uchar')');
	header.nSteps       = fread(fid,1,'int32');
	header.dt           = fread(fid,1,'float32');
	header.skip         = fread(fid,1,'int32');
	header.flagSv       = fread(fid,1,'int32');
	header.flagDm       = fread(fid,1,'int32');
     if(modelflag==0)
         header.zDes         = fread(fid,1,'float32');
     end;

if(modelflag == 0)
% function [header]=readVgrid_elcirc(fid,header,modelflag)
  	header.vgrid.startPos  = ftell(fid);
  	header.vgrid.zMsl      = fread(fid,1,'float32');
  	header.vgrid.nLevels   = fread(fid,1,'int32');
  	header.vgrid.zLevel    = fread(fid,header.vgrid.nLevels,'float32');
elseif(modelflag == 1)
%  function  [header]=readVgrid_selfe(fid,headermodelflag)
    header.vgrid.startPos = ftell(fid);
    header.vgrid.ivcor = fread(fid,1,'float32');
    header.vgrid.h_0 = fread(fid,1,'float32');
    header.vgrid.h_c = fread(fid,1,'float32');
    header.vgrid.theta_b = fread(fid,1,'float32');
    header.vgrid.theta_f = fread(fid,1,'float32');
    header.vgrid.zLevel = fread(fid,1,'int32');
    header.vgrid.sigma = fread(fid,[1 header.vgrid.zLevel],'float32');
end;
function [header]=readHgrid(fid,header,modelflag)
	header.hgrid.startPos   = ftell(fid);
    header.hgrid.type       = 'gr'; %it is a grid without a boundary
    header.hgrid.np         = fread(fid,1,'int32')
	header.hgrid.ne         = fread(fid,1,'int32')
		
	pos1=ftell(fid);
	hgridTmp                = fread(fid,[3 header.hgrid.np],'float32')';
	pos2=ftell(fid);
    
	header.hgrid.x          = hgridTmp(:,1);
	header.hgrid.y          = hgridTmp(:,2);
    header.hgrid.depth      = hgridTmp(:,3);
    if(modelflag == 0)
        header.hgrid.nodes      = [[1:header.hgrid.np]' hgridTmp(:,4)];
        fseek(fid, -(pos2-pos1), 'cof');
        hgridTmp                = fread(fid,[4 header.hgrid.np],'int32')';
        header.hgrid.bottomLayer = hgridTmp(:,4);
        if length(findstr(header.dataFormat,'v2.0'))
            elem                = fread(fid,header.grid.ne*3,'int32');
            header.hgrid.elem   = [[1:size(elem,1)]' elem];
        elseif length(findstr(header.dataFormat,'v3.0'))
            header=readElem_3(fid,header);
        else
            disp('couldnt match the data format to read elements')
        end
    elseif(modelflag == 1)
        header.hgrid.nodes = [[1:header.hgrid.np]' hgridTmp(:,3)];
        header=readElem_3(fid,header);
    end;


function [header]=readElem_3(fid,header,modelflag)
    ne                  = header.hgrid.ne;
    elem                = nan*ones(ne,5);
    for i=1:ne
        toq = fread(fid,1,'int32');
        elem(i,1)       = toq;
        if toq == 3
            elem(i,2:4) = fread(fid,3,'int32')';
        elseif toq == 4
            elem(i,2:5) = fread(fid,4,'int32')';
        else
            error(['oops problems with reading element ' num2str(i),' ',num2str(toq)])
        end
    end
    header.hgrid.elem   = [[1:size(elem,1)]' elem];;
    
    