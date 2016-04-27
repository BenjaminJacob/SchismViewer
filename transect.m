


h=imline


h = imfreehand



[xi yi]=ginput(2)


hold on
plot(xi,yi,'r--')



for i=1:length(xi)
    [was wo(i)]=min( (ll.hgrid.x-xi(i)).^2+(ll.hgrid.y-yi(i)).^2);
end


x1=gr.hgrid.x(wo(1));
y1=gr.hgrid.y(wo(1));


plot(gr.hgrid.x(wo),gr.hgrid.y(wo),'r--')


dx=diff(gr.hgrid.x(wo));
dy=diff(gr.hgrid.y(wo));
dy_dx=dy/dx;

L=sqrt(dx.^2+dy.^2);

ds=1000;

n=1+floor(L/ds);


dxx=dx/n;
dyy=dy/n;
xs=x1+(0:n)*dxx;
ys=y1+(0:n)*dyy;

p=[xs' ys'];
 

 
 gr_plot(gr.hgrid,gr.hgrid.depth)
 hold on
 plot(xs,ys,'ro')
 
 parents=find_parent(gr,p);

 
  