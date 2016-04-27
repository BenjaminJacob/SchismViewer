%generation of grid structures
clear
addpath '/home/users/frolovs/melio/'
gr=gr_readGrid([]);		%initialise an empty grid structure
gr=gr_readGrid(gr,['.'],'dr')	%read all grid files in the curent directory
gr=gr_tri(gr);			%triangualte files for visualisation
save gr gr			%save

