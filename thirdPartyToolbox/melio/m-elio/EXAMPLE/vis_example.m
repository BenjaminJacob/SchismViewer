%example of visualisation of salinity and hotstart files
%assumes that grid info is pregenerated in grid_db11

clear
addpath '/home/users/frolovs/melio/'
load grid_db11/gr

%%%%%%%%%%%%%%%%%%%%%%%%%
%read and vis salinity file
%read file
salt.h	= eb_readHeader( '1_salt.63' )	%read header of the salinity file
					%read first two timesteps of the salinity file
[salt.h,salt.ts,salt.data] = eb_readTimeStep(salt.h,[1 2]);
%convert elcirc binary data layout to layout similar to hotstart files, hotstart files are esear to visualise		
vis.h 	= salt.h; vis.data = salt.data(:,1); %visualise only the first timestep
vis.hts = map_eb2hts(vis)';

%plot layer 52 of salinity defined on the hgrid 
gr_plot(gr.hgrid,vis.hts(:,52),[0 34])

%%%%%%%%%%%%%%%%%%%%%%%%%%
%read and vis hotstart file
clear salt vis
%read hotstart ignore turbulence parameters
[hts]=hotstart_read_elcirc('hotstart.in',gr,1)
%visualize salinity defined at grid points
gr_plot(gr.hgrid,hts.snd(52,:)',[0 34])
%visualise salinity defined at side centers
gr_plot(gr.sideCent,hts.ssd(52,:)',[0 34])
%note the last figure looks like garbage because matlab's triangualator doesn't treat bounders correctly

