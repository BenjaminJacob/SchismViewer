function simplexid=intriangle(hgrid,p)



simplexid=nan(size(p,1),1);

i_tri=find(hgrid.elem(:,2)==3);

xs=hgrid.x;
ys=hgrid.y;


for ip=1:size(p,1)
    
    %split triangles and quads
    
    verts=hgrid.elem(i_tri,3:end-1);
    x2=xs(verts);
    y2=ys(verts);

    
    
    %gibt es punkte in 3 versch quadranten oder 2 in einem plus ein im ggü
    %liegenden                            %punkt in
    iq1=(x2-p(ip,1) >0 & y2-p(ip,2) >=0); %quadrant 1
    iq2=(x2-p(ip,1) <=0 & y2-p(ip,2) >0); %quadrant 2
    iq3=(x2-p(ip,1) < 0 & y2-p(ip,2)<=0); %quadrant 3
    iq4=(x2-p(ip,1) >= 0 & y2-p(ip,2)<0);  %quadrant 4
    
    %vorauswahl moegliche candidaten

    
    %% find intriangles
    i_tri=find(hgrid.elem(:,2)==3);
    verts=hgrid.elem(i_tri,3:end-1);
    x2=xs(verts);
    y2=ys(verts);

    
    ipos=find(...
        min(sum(iq1,2),1+logical(sum(iq3,2)))+min(sum(iq3,2),1+logical(sum(iq1,2))) + ...  %diagonale q1 und q3
        min(sum(iq2,2),1+logical(sum(iq4,2)))+min(sum(iq4,2),1+logical(sum(iq2,2)))...    %diagonale q2 un q4
        ==3);
    
    
    for k=1:length(ipos)
        in = inpolygon(p(ip,1),p(ip,2),x2(ipos(k),:),y2(ipos(k),:));
    
        %plot(xs(verts(ipos(k),:)),ys(verts(ipos(k),:)),'r')
        
        if in

            simplexid(ip)=i_tri(ipos(k));
            
            break
       
        end
    end
    
    %% find in quads
    
        i_quad=find(hgrid.elem(:,2)==4);
    verts=hgrid.elem(i_quad,3:end);
    x2=xs(verts);
    y2=ys(verts);

    
%     ipos=find(...
%         min(sum(iq1,2),1+logical(sum(iq3,2)))+min(sum(iq3,2),1+logical(sum(iq1,2))) + ...  %diagonale q1 und q3
%         min(sum(iq2,2),1+logical(sum(iq4,2)))+min(sum(iq4,2),1+logical(sum(iq2,2)))...    %diagonale q2 un q4
%         ==4);
    
    
    for k=1:length(i_quad)
        in = inpolygon(p(ip,1),p(ip,2),x2(k,:),y2(k,:));
    
        %plot(xs(verts(ipos(k),:)),ys(verts(ipos(k),:)),'r')
        
        if in
            simplexid(ip)=i_quad(k);
            
            break
       
        end
    end

end