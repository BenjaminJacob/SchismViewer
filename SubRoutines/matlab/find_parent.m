
function parents=find_parent(gr,p)

parents=zeros(size(p,1),1);

for ip=1:size(p,1)
    
    [dists isorted]=sort(sqrt( (gr.hgrid.x-p(ip,1)).^2+(gr.hgrid.y-p(ip,2)).^2));
    
    found=0;
    for i=1:20
        
        parents_pot=find(nansum(gr.hgrid.elem(:,3:6)==isorted(i) ,2));
        
        
        for parent=parents_pot'
            
            if gr.hgrid.elem(parent,2)==3 %tris
                xv=gr.hgrid.x(gr.hgrid.elem(parent,[3:5 3]));
                yv=gr.hgrid.y(gr.hgrid.elem(parent,[3:5 3]));
            else %quad
                xv=gr.hgrid.x(gr.hgrid.elem(parent,[3:6 3]));
                yv=gr.hgrid.y(gr.hgrid.elem(parent,[3:6 3]));
            end
            
            
            
            inside=inpolygon(p(ip,1),p(ip,2),xv,yv);
            
            if inside
            
                found=1;
                parents(ip)=parent;
                break
            end

        end
        
        if found
           break 
        end
        
    end
end

