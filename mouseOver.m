
function woop
figure()
set (gcf, 'WindowButtonMotionFcn', @mouseover);

function [data] = mouseover(gcbo,eventdata,handles)
keyboard    
c = get (gca, 'CurrentPoint'); % get mouse coordinates
    
