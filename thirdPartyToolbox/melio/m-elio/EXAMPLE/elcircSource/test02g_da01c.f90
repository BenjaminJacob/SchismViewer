!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                	!
!                                ELCIRC:                                    		!
!         A Three-Dimensional Baroclinic Model for Unstructured Grids         		!
!		        Version 5.01.01c (Sept. 05, 2003)            	      		!
!                                                                             		!
!                 Center for Coastal and Land-Margin Research                 		!
!             Department of Environmental Science and Engineering             		!
!                   OGI School of Science and Engineering,		      		!
!	              Oregon Health & Science University		      		!
!                       Beaverton, Oregon 97006, USA                            	!
!								 	      		!
!                   Scientific direction: Antonio Baptista                    		!
!                   Code development: Joseph Zhang (since April 2001).	      		!
!		                      Ed Myers (1999-May 2001);               		!
!											!
!	        Copyright 1999-2003 Oregon Health and Science University		!
!  		               All Rights Reserved					!
!                                                                             		!
!     	Formulation of Elcirc was inspired by Casulli and Zanolli (1998).       	!
!	Sparse matrix solver dsrc2c.f90 was adapted from ITPACK:			!
!	D. Young and D. Kincaid. ``The ITPACK Package for Large Sparse Linear Systems,''!
!	in Elliptic Problem Solvers, (M. Schultz, ed.), Academic Press, 		!
!	New York, 1981, pp.163-185. 							!
!    									      		!
! 	The heat exchange module makes use of the bulk aerodynamic surface flux 	!
!	algorithm introduced by Zeng et al (1998), and the polynomial fits to 		!
!	saturation vapor pressure of Flatau et al (1992):				! 
!	Zeng, X., M. Zhao, and R. E. Dickinson, 1998:  Intercomparison of bulk		!
!	aerodynamic algorithms for the computation of sea surface fluxes using		!
!	TOGA COARE and TAO data.  J. Clim., 11, 2628-2644.				!
!	Flatau, P. J., R. L. Walko and W. R. Cotton, 1992:  Polynomial fits to		!
!	saturation vapor pressure.  J. Appl. Meteor., 31, 1507-1513.			!
!											!
!	Attenuation of solar radiation (and solar heating) within the water column	!
!	is based upon the expression given by Paulson and Simpson (1977), for the	!
!	water types defined by Jerlov (1968):						!
!	Jerlov, N. G., Optical Oceanography, Elsevier, 1968.				!
!	Paulson, C. A., and J. J. Simpson, Irradiance measurements in the upper		!
!	ocean, J. Phys. Oceanogr., 7, 952-956, 1977.					!
!											!
!	In addition, the module must be linked with Hierarchical Data Format (HDF)	!
!	libraries developed at the National Center for Supercomputing Applications	!
!	(NCSA) at University of Illinois at Urbana-Champaign (UIUC).  See		!
!	http://hdf.ncsa.uiuc.edu/ for additonal information, source code, etc.		!
!											!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             		!

!...  Data type consts
      module kind_par
        implicit none
        integer, parameter :: sng_kind1=4
        integer, parameter :: dbl_kind1=8
      	real(kind=dbl_kind1), parameter :: small1=1.e-6 !small non-negative number; must be identical to that in global
      end module kind_par


!...  definition of variables
!...
!
!************************************************************************
!     			mnp < mne < mns					*
!************************************************************************
!
      module global
      	implicit none
	integer, parameter :: sng_kind=4
        integer, parameter :: dbl_kind=8

!...  	Dimensioning parameters
        integer, parameter :: mnp=35000
      	integer, parameter :: mne=51000
      	integer, parameter :: mns=85000
      	integer, parameter :: mnv=90
      	integer, parameter :: mnei=30 !neighbor
      	integer, parameter :: mnope=6 !# of open bnd segements
       	integer, parameter :: mnond=20000 !max. # of open-bnd nodes on each segment
      	integer, parameter :: mnland=50 !# of land bnd segements
      	integer, parameter :: mnlnd=10000 !max. # of land nodes on each segment
      	integer, parameter :: mnoe=20000 !max. # of open-bnd elements on each segment
      	integer, parameter :: mnosd=20000 !max # of open-bnd sides on each segment
      	integer, parameter :: mnbfr=15 !# of forcing freqs.
      	integer, parameter :: itmax=5000 !# of iteration for itpack solvers used for dimensioning
      	integer, parameter :: nwksp=6*mne+4*itmax !available work space for itpack solvers
      	integer, parameter :: nbyte=4
      	integer, parameter :: mnout=100 !max. # of output files
      	integer, parameter :: mirec=1109000000 !997000000) !max. record # to prevent output ~> 4GB
      	real(kind=dbl_kind), parameter :: small1=1.e-6 !small non-negative number; must be identical to that in kind_par

!...  	Important variables
      	integer :: np,ne,ns,nvrt,ieqstate
      	real(kind=dbl_kind) :: zmsl,h0,denom,q2min,rho0,dt
!	Consts. used in UB closure
        character(len=2) :: mid,stab
        real(kind=dbl_kind) :: ubd0,ubd1,ubd2,ubd3,ubd4,ubd5,ubs0,ubs1,ubs2,ubs4,ubs5,ubs6, &
     &a2_cm03,schk,schpsi

!...    Output handles
        character(len=48) :: start_time,version,data_format='DataFormat v3.0'
        character(len=12) :: ifile_char
        character(len=48), dimension(mnout) :: outfile,variable_nm,variable_dim
        integer :: ihot,nrec,nspool,igmp,noutgm,ifile,noutput
        integer, dimension(mnout) :: ichan,irec,iof
        real(kind=dbl_kind), dimension(mnout) :: vpos
        
! evm   character buffers for binary type conversions
	integer :: iwrite
        character(len=48)   :: a_48
        character(len=16)   :: a_16
        character(len=8)    :: a_8
        character(len=4)    :: a_4


!...  1D arrays
	integer :: nne(mnp),nnp(mnp),kfs(mns),kbs(mns),kbp(mnp),kfp(mnp), &
	   &kfe(mne),kbe(mne)

	real(kind=dbl_kind), dimension(mnp) :: x,y,dp,peta
	real(kind=dbl_kind), dimension(mne) :: areas,radiel,xctr,yctr,dpe,eta1,eta2
	real(kind=dbl_kind), dimension(mns) :: snx,sny,distj,xcj,ycj,delj,dps,xlmin1,xlmin2

	real(kind=dbl_kind) :: delz(mnv),ztot(0:mnv),tem0(mnp,mnv),sal0(mnp,mnv)

!...  2D arrays
	integer :: i34(mne),nm(mne,4),nx(4,4,3),ic3(mne,4),ine(mnp,mnei),js(mne,4),is(mns,2), &
	   &isidenode(mns,2),inp(mnp,mnei),iself(mnp,mnei)

	real(kind=dbl_kind) :: ssign(mne,4),dz(mns,mnv),dzhalf(mns,0:mnv),dc(mne,mnv), &
	   &dchalf(mne,mnv),dzp(mnp,mnv),dzphalf(mnp,mnv),we(mne,0:mnv), &
	   &vn2(mns,0:mnv),vt2(mns,0:mnv),tsd(mns,mnv),ssd(mns,mnv),tnd(mnp,mnv),snd(mnp,mnv), &
	   &srho(mns,mnv),erho(mne,mnv),prho(mnp,mnv),q2(mns,0:mnv),xl(mns,0:mnv)

	real(kind=dbl_kind), dimension(mnp,0:mnv) :: uu1,vv1,ww1, uu2,vv2,ww2

!...  serg
!...  declaration of da varaibles in mod global
        character(len=48), dimension(mnout) :: da_outfile,da_variable_nm,da_variable_dim
        integer :: da_noutput
        integer, dimension(mnout) :: da_ichan, da_iof, da_irec
        real(kind=dbl_kind), dimension(mnout) :: da_vpos
!... other usefull vars
        real(kind=dbl_kind) :: da_ugsc(mns,mnv), da_vgsc(mns,mnv)
        integer :: da_ihot_type
!... serg 

      end module global

!...  Main program
      program elcirc
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

!...  Output handles
      real(kind=sng_kind) :: floatout,floatout2,st,en
      character(len=12) :: it_char
      character(len=40) :: date,timestamp

!... Geometry
      dimension xlon(mnp),ylat(mnp),x1j(mns),y1j(mns),x2j(mns),y2j(mns)
      dimension xlon_e(mne),ylat_e(mne)

!...  Boundary forcings
      character(len=5) :: bountag(mnbfr)
      character(len=10) :: fq_nm(mnbfr)
      dimension isbnd(mnp),sarea(mnp),ibad(mnp),idrynode(mnp)
      dimension isbnode(mnp,2),isbe(mne),isbs(mns)
      dimension nond(mnope),iond(mnope,mnond)
      dimension nlnd(mnland),ilnd(mnland,mnlnd)
      dimension nosd(mnope),iosd(mnope,mnosd),noe(mnope),ioe(mnope,mnoe)
      dimension isflowside(mns),isflowside2(mns),cwidth(mnope)
      dimension iettype(mnope),ifltype(mnope)
      dimension itetype(mnope),isatype(mnope)
      dimension tamp(mnbfr),tnf(mnbfr),tfreq(mnbfr),jspc(mnbfr),tear(mnbfr)
      dimension fun_lat(mnp,0:2),fun_lat_e(mne,0:2)
      dimension amig(mnbfr),ff(mnbfr),face(mnbfr)
      dimension emo(mnope,mnbfr,mnoe),efa(mnope,mnbfr,mnoe)
      dimension vmo(mnope,mnbfr,mnosd),vfa(mnope,mnbfr,mnosd)
      dimension eth(mnope),qth(mnope),tth(mnope),sth(mnope)
      dimension uth(mnope),vth(mnope),ath(mnope)

!...  Flow arrays
      dimension tsdbt(mns,mnv),ssdbt(mns,mnv),vnbt(mns,mnv)
      dimension tndbt(mnp,mnv),sndbt(mnp,mnv),vtbt(mns,mnv)
      dimension zhat1(mns,mnv),ghat(mns,mnv),gvec(mnv)
      dimension zhat2(mns,mnv),fhat(mns,mnv),fvec(mnv),ndelt(mnp,mnv)
      dimension icoef(mne+1),jcoef(mne*5),e2coef(mne*5),qel(mne)
      dimension sparsem(mne,0:4),elbc(mne),imape(mne),qel2(mne),eta3(mne)
      dimension windx1(mnp),windy1(mnp),windx2(mnp),windy2(mnp)
      dimension windx(mns),windy(mns),tau_n(mns),tau_t(mns),advt(mns)
      dimension pr1(mnp),airt1(mnp),shum1(mnp),pr2(mnp),airt2(mnp),shum2(mnp)
      dimension sflux(mnp),srad(mnp)
      dimension fluxsu(mnp),fluxlu(mnp),hradu(mnp),hradd(mnp)
      dimension tauxz(mnp),tauyz(mnp),windxp(mnp),windyp(mnp)
      dimension tk(mns),cori(mns),bdragc(mns),Cd(mns),roughness(mns)
      dimension vdiff(mns,mnv),tdiff(mns,mnv)
      dimension tdiffp(mnp,mnv),etaic(mne),relax(mnp)
      dimension rich(mns,mnv)
      dimension hor_t(mne),hor_s(mne)

!...  Wild-card arrays
      dimension nwild(mnp),nwild2(mne),swild(mnp) !,slope(mne)

!... Solver arrays for TRIDAG
      dimension alow(mnv),bdia(mnv),cupp(mnv),rrhs(mnv,5)
      dimension soln(mnv,5),gam(mnv)
      dimension trhs(mns,mnv),srhs(mns,mnv)

!     MY-G turbulence closure arrays
      dimension qdiff(mns,mnv),qdiff2(mns,mnv),q2tmp(mnv),xltmp(mnv)
      dimension q2bt(mns,mnv),xlbt(mns,mnv)
      dimension rzbt(mns,mnv),shearbt(mns,mnv),xlmax(mnv),diffmax(mns)
      dimension q2nd(mnp,mnv),xlnd(mnp,mnv)
!     Test arrays
      dimension testa(mns)

!...  variables used by the itpack solvers
      dimension iwksp(3*mne),wksp(nwksp),iparm(12),rparm(12)


!...
!...  First executible statement of Elcirc
!...

!... Initialize arrays
      idrynode=0 !not dry
      isbnd=0
      isbe=0
      isbs=0
      pr1=0; pr2=0 !uniform pressure (the const. is unimportant)
!	for output
      airt1=0; shum1=0;  airt2=0; shum2=0; srad=0; fluxsu=0; fluxlu=0 
      hradu=0; hradd=0; sflux=0; windxp=0; windyp=0

!...  define some constants and initial values
!...
      omega=7.29d-5 !angular freq. of earth rotation 
      rearth=6378206.4 !earth radius
      g=9.81
      rho0=1000. !1025. !ref. density for S=33 and T=10C
      pi=3.141592653589793d0
      deg2rad=pi/180
      rad2deg=180/pi
      shw=4184  !specific heat of pure water

      do k=3,4
        do i=1,k
          do j=1,k-1
	    nx(k,i,j)=i+j
      	    if(nx(k,i,j)>k) nx(k,i,j)=nx(k,i,j)-k
	    if(nx(k,i,j)<1.or.nx(k,i,j)>k) then
	      write(*,*)'nx wrong',i,j,k,nx(k,i,j)
	      stop
	    endif
          enddo !j
        enddo !i
      enddo !k


!...  Boundary layer characteristics
      hestu=60. !threshold depth for estuary (m)
      hcont=100.  !threshold depth for continuental shelf
      bthick_estu=0.25 !char. BL thickness for estuary
      bthick_cont=1.00 !char. BL thickness for continuental shelf
      bthick_ocea=29.8 !char. BL thickness for open ocean
      if(hestu.gt.hcont) then
        write(*,*)'Check depth scales',hestu,hcont
      stop
      endif

!...  set the maximum number of digits of precision the grid can be expected to have
!...  note: if the grid was build on a 32 bit computer, it should be accurate to about
!...  7 digits.  Thus, if the nodal spacing requires more than 5 digits of precision,
!...  it is unlikely that the model results will be trustworthy.

      nprec=5



!                                                                             *
!******************************************************************************
!                                                                             *
!			open input files				      *
!                                                                             *
!******************************************************************************
!                                                                             *

      open(14,file='hgrid.gr3',status='old')
      open(15,file='param.in',status='old')
      open(19,file='vgrid.in',status='old')

      open(11,file='fort.11') !fatal error message output
      open(12,file='fort.12') !non-fatal error message output
      open(16,file='mirror.out')
      rewind(11)
      rewind(16)

      call date_and_time(date,timestamp)
      write(16,*)'Run begins at ',date,timestamp


!...  read the vertical layers information from vgrid.in
!...
      read(19,*) nvrt,zmsl
      if(nvrt.gt.mnv) then
 	write(11,*)'nvrt > mnv'
      	stop
      endif

      ztot(0)=0 !0th level
      delzmin=1.0d15 !min. delz(i)
      do i=1,nvrt
        read(19,*) jki,delz(i),ztot(i)
      	if(delz(i).lt.delzmin) delzmin=delz(i)
      	if(dabs(ztot(i)-zmsl).lt.1.0d-5) then
          write(11,*)'zmsl is too close to level line',i
          stop
      	endif
      enddo
      close(19)

!     Check delz and ztot
      do i=1,nvrt
      	if(dabs(ztot(i-1)+delz(i)-ztot(i)).gt.1.e-6) then
          write(11,*)'Wrong input fort.19 at level',i,ztot(i-1),delz(i)
          stop
     	endif
      enddo
      if(zmsl.gt.ztot(nvrt).or.zmsl.lt.0) then
      	write(11,*)'MSL above/below  all levels'
      	stop
      endif

!...  read unit 15 file
!...

      read(15,'(a48)') version
      read(15,'(a48)') start_time
      read(15,*) ipre !pre-processing flag (to output obe.out)
      read(15,*) nscreen
      read(15,*) iwrite !write mode

      read(15,*) ihot
      if(ihot.lt.0.or.ihot.gt.2) then
      	write(11,*)'Unknown ihot',ihot
      	stop
      endif
!... serg
!... type of the hotstart file elcirc or da
      read(15,*) da_ihot_type
      if(da_ihot_type.ne.1.and.da_ihot_type.ne.2) then
        write(11,*)'Unknown da_ihot_type',da_ihot_type
        stop
      endif    
!... serg

      read(15,*) ics
      if(ics.ne.1.and.ics.ne.2) then
        write(11,*)'Unknown ics',ics
      	stop
      endif

!...  Center of projection in degrees
      read(15,*) slam0,sfea0
      slam0=slam0*deg2rad
      sfea0=sfea0*deg2rad

!... Enough info to read unit 14 grid file
!...
      read(14,*) 
      read(14,*) ne,np
      if(ne.gt.mne.or.np.gt.mnp) then
        write(11,*)'Increase mne/mnp',mne,mnp,ne,np
        stop
      endif

      do i=1,np
        if(ics.eq.1) then
          read(14,*) j,x(i),y(i),dp(i)
        else if(ics.eq.2) then
          read(14,*) j,xlon(i),ylat(i),dp(i)
          ylat(i)=ylat(i)*deg2rad
          xlon(i)=xlon(i)*deg2rad
          call cpp(x(i),y(i),xlon(i),ylat(i),slam0,sfea0)
        endif
	if(dp(i).le.0) idrynode(i)=1
      enddo !i=1,np

      do i=1,ne
        read(14,*) j,i34(i),(nm(i,k),k=1,i34(i))
	if(i34(i)/=3.and.i34(i)/=4) then
	  write(*,*)'Unknown type of element',i
	  stop
	endif

	n1=nm(i,1)
	n2=nm(i,2)
	n3=nm(i,3)
	if(i34(i)==3) then
          areas(i)=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
	else !quad
	  n4=nm(i,4)
!         Check convexity
          ar1=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
          ar2=signa(x(n1),x(n3),x(n4),y(n1),y(n3),y(n4))
          ar3=signa(x(n1),x(n2),x(n4),y(n1),y(n2),y(n4))
          ar4=signa(x(n2),x(n3),x(n4),y(n2),y(n3),y(n4))
          if(ar1<=0.or.ar2<=0.or.ar3<=0.or.ar4<=0) then
            write(*,*)'Concave quadrangle',i,ar1,ar2,ar3,ar4
            stop
          endif

          areas(i)=ar1+ar2
	endif
        if(areas(i).le.0.0) then
          write(11,*)'Negative area at',i
          stop
        endif
        radiel(i)=dsqrt(areas(i)/pi)  !equivalent radius
      enddo !i=1,ne

!     Open bnds
      read(14,*) nope
      if(nope.gt.mnope) then
      	write(11,*) 'nope > mnope' 
        stop
      endif

      read(14,*) neta
      jnmm=0
      do k=1,nope
        read(14,*) nond(k)
        if(nond(k).gt.mnond) then
  	  write(11,*) 'nond(k) > mnond'
          stop
        endif
        do i=1,nond(k)
          read(14,*) iond(k,i)
          isbnd(iond(k,i))=k
        enddo
        jnmm=jnmm+nond(k)
      enddo

      if(neta.ne.jnmm) then
        write(11,*)'neta /= total # of open bnd nodes',neta,jnmm
        stop
      endif

!     Land bnds
      read(14,*) nland
      if(nland.gt.mnland) then
      	write(11,*) 'nland > mnland'
        stop
      endif

      read(14,*) nvel
      do k=1,nland
        read(14,*) nlnd(k)
        if(nlnd(k).gt.mnlnd) then
          write(11,*)'nlnd(k) > mnlnd',nlnd(k),mnlnd
          stop
        endif
        do i=1,nlnd(k)
          read(14,*) ilnd(k,i)
          if(isbnd(ilnd(k,i)).eq.0) isbnd(ilnd(k,i))=-1 !overlap of open bnd
        enddo
      enddo !k=1,nland
      close(14)
!...  End fort.14

!                                                                             *
!                                                                             *
!******************************************************************************
!                                                                             *
!     			Compute geometry 				      *
!                                                                             *
!******************************************************************************
!                                                                             *
!                                                                             *

!...  compute the elements connected to each node 
      do i=1,np
        nne(i)=0
      	nnp(i)=0
      	sarea(i)=0
      enddo

      do i=1,ne
        do j=1,i34(i)
          nd=nm(i,j)
          nne(nd)=nne(nd)+1
          if(nne(nd).gt.mnei) then
            write(11,*)'Too many neighbors',nd
            stop
          endif
          ine(nd,nne(nd))=i
          iself(nd,nne(nd))=j
          sarea(nd)=sarea(nd)+areas(i)
        enddo
      enddo

!...  Check hanging nodes
!...
      ihang=0
      do i=1,np
        if(nne(i).eq.0) then
          ihang=1
          write(11,*)'Hanging node',i
        endif
      enddo
      if(ihang.eq.1) then
        write(11,*)'Check fort.11 for hanging nodes'
        stop
      endif
      
!...  Compute surrounding nodes
!...
      do i=1,np
      	do j=1,nne(i)
          ie=ine(i,j)
          loopk: do k=1,i34(ie)
            nd=nm(ie,k)
            do l=1,nnp(i)
              if(nd.eq.i.or.nd.eq.inp(i,l)) cycle loopk
            enddo !l
            nnp(i)=nnp(i)+1
            if (nnp(i).gt.mnei) then
              write(11,*)'Too many surrounding nodes',i
              stop
            endif
            inp(i,nnp(i))=nd
	  end do loopk !k
      	enddo !j
      enddo !i=1,np

!     Compute ball info
      do i=1,ne
        do j=1,i34(i)
          ic3(i,j)=0 !index for bnd sides
          nd1=nm(i,nx(i34(i),j,1))
          nd2=nm(i,nx(i34(i),j,2))
          do k=1,nne(nd1)
            iel=ine(nd1,k)
            if(iel/=i.and.(nm(iel,1)==nd2.or.nm(iel,2)==nd2.or.&
     &nm(iel,3)==nd2.or.(i34(iel)==4.and.nm(iel,4)==nd2))) ic3(i,j)=iel
          enddo !k
        enddo !j
      enddo !i

!...  compute the sides information
!...
      ns=0 !# of sides
      do i=1,ne
        do j=1,i34(i)
          nd1=nm(i,nx(i34(i),j,1))
          nd2=nm(i,nx(i34(i),j,2))
          if(ic3(i,j)==0.or.i<ic3(i,j)) then !new sides
            ns=ns+1 
            if(ns.gt.mns) then
              write(11,*)'Too many sides'
              stop
            endif
            js(i,j)=ns
            is(ns,1)=i
            isidenode(ns,1)=nd1
            isidenode(ns,2)=nd2
            x1j(ns)=x(nd1)
            y1j(ns)=y(nd1)
            x2j(ns)=x(nd2)
            y2j(ns)=y(nd2)
            xcj(ns)=(x1j(ns)+x2j(ns))/2
            ycj(ns)=(y1j(ns)+y2j(ns))/2
            dps(ns)=(dp(nd1)+dp(nd2))/2
            distj(ns)=dsqrt((x2j(ns)-x1j(ns))**2+(y2j(ns)-y1j(ns))**2)
            if(distj(ns).eq.0) then
              write(11,*)'Zero side',ns
              stop
            endif
            is(ns,2)=ic3(i,j) !bnd element => bnd side
!	Corresponding side in element ic3(i,j)
            if(ic3(i,j)/=0) then !old internal side
              iel=ic3(i,j)
              index=0
              do k=1,i34(iel)
                if(ic3(iel,k)==i) then
                  index=k
                  exit
                endif
              enddo !k
   	      if(index==0) then
                write(*,*)'Wrong ball info',i,j
                stop
              endif
              js(iel,index)=ns
            endif !ic3(i,j).ne.0
          endif !ic3(i,j)==0.or.i<ic3(i,j)
        enddo !j=1,i34(i)
      enddo !i=1,ne

      if(ns.lt.ne.or.ns.lt.np) then
        write(11,*)'Weird grid with ns < ne or ns < np',np,ne,ns
        stop
      endif

      do i=1,ne
        do j=1,i34(i)
          jsj=js(i,j)
          ssign(i,j)=(is(jsj,2)-2*i+is(jsj,1))/(is(jsj,2)-is(jsj,1))
        enddo
      enddo

      if(nscreen.eq.1) write(*,*) 'There are',ns,' sides in the grid...'
      write(16,*) 'There are',ns,' sides in the grid...'

!...  Compute open bnd elements and sides
!...  Each bnd node is assigned 2 flags to account for overlap of 2
!...  open bnds.
      do i=1,np
        isbnode(i,1)=0
        isbnode(i,2)=0
      enddo
      do i=1,nope
        do j=1,nond(i)
          nd=iond(i,j)
          if(isbnode(nd,1).eq.0) then
            isbnode(nd,1)=i
          else if(isbnode(nd,2).eq.0) then
            isbnode(nd,2)=i
          else
            write(11,*)'Illegal open bnd at node',nd
            stop
          endif
        enddo !j
      enddo !i=1,nope

      do j=1,nope
        nosd(j)=0 !initialize
        noe(j)=0
      enddo

      loop1: do i=1,ns
      	if(is(i,2)==0) then
          nd1=isidenode(i,1)
          nd2=isidenode(i,2)
          do j=1,nope
            if((isbnode(nd1,1).eq.j.or.isbnode(nd1,2).eq.j).and.&
     &         (isbnode(nd2,1).eq.j.or.isbnode(nd2,2).eq.j)) then
              nosd(j)=nosd(j)+1
              if(nosd(j).gt.mnosd) then
                write(11,*)'Too many open sides',j
                stop
              endif
              iosd(j,nosd(j))=i
              cycle loop1 !i can only be in one open bnd segment
            endif
          enddo !j=1,nope
        endif !is(i,2)=0
      end do loop1 !i=1,ns

      loop2: do i=1,ne
        do j=1,i34(i)
          nd=nm(i,j)
          do k=1,nope
            if(isbnode(nd,1)==k.or.isbnode(nd,2)==k) then
              noe(k)=noe(k)+1
              if(noe(k).gt.mnoe) then
                write(11,*)'Too many open elements',k
                stop
              endif
              ioe(k,noe(k))=i
!	In the case of dual role, choose any one of the segment
              cycle loop2
            endif
          enddo !k
        enddo !j
      end do loop2 !i=1,ne

!...  compute centers of each element & dpe
!...
      do i=1,ne
	xctr(i)=0
	yctr(i)=0
	dpe(i)=-5.e8
	do j=1,i34(i)
          xctr(i)=xctr(i)+x(nm(i,j))/i34(i)
          yctr(i)=yctr(i)+y(nm(i,j))/i34(i)
	  if(dps(js(i,j))>dpe(i)) dpe(i)=dps(js(i,j))
	enddo !j
      enddo !i=1,ne

!...  Output obe.out, centers.bp and sidecenters.bp if ipre.ne.0
      if(ipre.ne.0) then
        open(32,file='obe.out')
        if(iwrite.eq.0) then
          write(32,*) nope,' # of open bnd'
          write(32,*) 'Element list:'
        else !evm
          write(32,"(i12,a)",advance="no") nope,' # of open bnd\n'
          write(32,"(a)",advance="no") 'Element list:\n'
        endif
        do i=1,nope
          if(iwrite.eq.0) then
            write(32,*) noe(i),' bnd #',i 
          else !evm
            write(32,"(i12,a,i12,a)",advance="no") noe(i),' bnd #',i,'\n'
          endif
          do j=1,noe(i)
            if(iwrite.eq.0) then
              write(32,*) j,ioe(i,j)
            else !evm
              write(32,"(2i12,a)",advance="no") j,ioe(i,j),'\n'
            endif
          enddo !j
        enddo !i
        if(iwrite.eq.0) then
          write(32,*) 'Side list:'
        else !evm
          write(32,"(a)",advance="no") 'Side list:\n'
        endif
        do i=1,nope
          if(iwrite.eq.0) then
            write(32,*) nosd(i),' bnd #',i 
          else !evm
            write(32,"(i12,a,i12,a)",advance="no") nosd(i),' bnd #',i,'\n'
          endif
          do j=1,nosd(i)
            if(iwrite.eq.0) then
              write(32,*) j,iosd(i,j)
            else !evm
              write(32,"(2i12,a)",advance="no") j,iosd(i,j),'\n'
            endif
          enddo
        enddo !i
        close(32)

        open(32,file='sidecenters.bp')
        if(iwrite.eq.0) then
          write(32,*) 'Sidegrid'
          write(32,*) ns
        else !evm
          write(32,"(a)",advance="no") 'Sidegrid\n'
          write(32,"(i12,a)",advance="no") ns,'\n'
        endif
        do i=1,ns
          if(iwrite.eq.0) then
            write(32,*) i,xcj(i),ycj(i),real(dps(i))
          else !evm
            write(32,"(i12,a,f19.9,a,f19.9,a,f12.6,a)",advance="no") i, &
     &" ",xcj(i),"      ",ycj(i),"    ",real(dps(i)),'\n'
          endif
        enddo !i
        close(32)

!... serg new output in preproccessing
!... ouptut connectivity for sidecenters
        open(32,file='sidecenters_conn.bp')
        if(iwrite.eq.0) then
          write(32,*) 'Sidegrid conn'
          write(32,*) ns
        else !evm
          write(32,"(a)",advance="no") 'Sidegrid\n'
          write(32,"(i12,a)",advance="no") ns,'\n'
        endif
        do i=1,ns
          if(iwrite.eq.0) then
             write(32,*) i, isidenode(i,1),isidenode(i,2), 1
          else !evm
            write(32,"(i12,a,i12,a,i12,a,i12,a)",advance="no") i, &
     &" ",isidenode(i,1),"      ",isidenode(i,2),"    ",1,'\n'
          endif
        enddo !i
        close(32)

!... serg end

        open(32,file='centers.bp')

        if(iwrite.eq.0) then
          write(32,*) 'centers pts'
          write(32,*) ne
        else !evm
          write(32,"(a)",advance="no") 'centers pts\n'
          write(32,"(i12,a)",advance="no") ne,'\n'
        endif
        do i=1,ne
          if(iwrite.eq.0) then
            write(32,*) i,xctr(i),yctr(i),real(dpe(i))
          else !evm
            write(32,"(i12,a,f19.9,a,f19.9,a,f12.6,a)",advance="no") i, &
     &" ",xctr(i),"      ",yctr(i),"    ",real(dpe(i)),'\n'
          endif
        enddo !i
        close(32)

        write(*,*)'obe.out successfully generated; bye'
        write(16,*)'obe.out successfully generated; bye'
        stop
      endif !ipre.ne.0

!...  compute angles (orientation of local x) and perpendicular distances for sides 
!...
      do j=1,ns
        if(is(j,2)/=0) then
          delj(j)=dsqrt((xctr(is(j,2))-xctr(is(j,1)))**2+(yctr(is(j,2))-yctr(is(j,1)))**2)
          if(delj(j).eq.0) then
            write(*,*)'Zero distance between centers',is(j,1),is(j,2)
            stop
          endif
        endif
        thetan=datan2(x1j(j)-x2j(j),y2j(j)-y1j(j))
      	snx(j)=dcos(thetan)
      	sny(j)=dsin(thetan)
      enddo

      if(nscreen.eq.1) write(*,*)'done computing geometry...'
      write(16,*)'done computing geometry...'

!...  Compute channel widths for type I b.c.
      do i=1,nope
        cwidth(i)=0
        do j=1,nosd(i)
          cwidth(i)=cwidth(i)+distj(iosd(i,j))
        enddo !j
      enddo !i

!...  classify boundary elements & sides
!...
      do i=1,nope
        do j=1,noe(i)
          isbe(ioe(i,j))=i
        enddo !j
        do j=1,nosd(i)
          isbs(iosd(i,j))=i
        enddo !j=1,nosd(i)
      enddo !i=1,nope

      do i=1,ns
      	if(is(i,2)==0.and.isbs(i)==0) isbs(i)=-1
      enddo

      if(nscreen.eq.1) write(*,*)'done classifying boundaries...'
      write(16,*)'done classifying boundaries...'

!...  Continue to read fort.15
!...

!... Implicitness for momentum
      read(15,*)thetai

!... Baroclinic flags
      read(15,*) ibc,ibtp
      if(ibc.ne.0.and.ibc.ne.1) then
      	write(*,*)'Unknown ibc'
      	stop
      endif
      if(ibtp.ne.0.and.ibtp.ne.1) then
      	write(*,*)'Unknown ibtp'
     	stop
      endif

      if(ibc.eq.0) then
        write(*,*)'You are using baroclinic model'
        write(16,*)'You are using baroclinic model'
      	read(15,*) nrampbc,drampbc 
      else !ibc=1
        if(ibtp.eq.0) then
          write(*,*)'Barotropic model without ST calculation'
          write(16,*)'Barotropic model without ST calculation'
        else !ibtp=1
          write(*,*)'Barotropic model with ST calculation'
          write(16,*)'Barotropic model with ST calculation'
        endif
      endif

      read(15,*) tempmin,tempmax,saltmin,saltmax
      if(tempmin.lt.0.or.tempmax.gt.40.or.saltmin.lt.0.or.saltmax.gt.42) then
  	write(*,*)'Specified ST range invalid'
      	stop
      endif

      read(15,*) rnday

!... dramp not used if nramp=0
      read(15,*) nramp,dramp
      if(nramp.ne.0.and.nramp.ne.1) then
      	write(*,*)'Unknown nramp',nramp
      	write(11,*)'Unknown nramp',nramp
        stop
      endif

      read(15,*) dt

!...  compute total number of time steps nt
      nt=rnday*86400.d0/dt+0.5

!...  input info on backtracking
!...
      read(15,*) nsubfl !flag
      if(nsubfl.lt.0.or.nsubfl.gt.2) then
        write(11,*)'Unknown subdivision flag',nsubfl
      	stop
      endif

      if(nsubfl.eq.0) then
        read(15,*) nsub
      	do i=1,mnp
          do j=1,nvrt
            ndelt(i,j)=nsub
          enddo
      	enddo
      else if(nsubfl.eq.1) then
      	open(41,file='ndelt.gr3',status='old')
        read(41,*)
        read(41,*) ntmp,npbp
        do i=1,npbp
          read(41,*)j,xtmp,ytmp,sub
          if(sub.le.0) then
            write(11,*)'# of subdivisions in btrack <=0',sub,i
            stop
          endif
          do j=1,nvrt
            ndelt(i,j)=sub
          enddo
        enddo !i
      	close(41)
      else !nsubfl=2
        read(15,*) ndelt_min,ndelt_max
      	if(ndelt_min.le.0.or.ndelt_min.gt.ndelt_max) then
          write(11,*)'Incorrect input for variable tracking!'
          write(11,*)ndelt_min,ndelt_max
          stop
        endif
      endif

!... Advection flag for momentum eq.
      read(15,*) nadv !flag
      if(nadv.ne.0.and.nadv.ne.1) then
        write(11,*)'Unknown advection flag',nadv
      	stop
      endif

      if(nadv.eq.0) then
      	open(42,file='adv.gr3',status='old')
        read(42,*)
        read(42,*) !ne,np
        do i=1,np
          read(42,*)j,xtmp,ytmp,tmp
          nwild(i)=tmp
      	enddo
        do i=1,ns
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          advt(i)=max0(nwild(n1),nwild(n2))
        enddo
        close(42)
      endif

!...  Minimum depth allowed
      read(15,*) h0 
      if(h0.le.0) then
        write(11,*)'h0 must be positive'
        stop
      endif
      denom=dmin1(h0*(0.5-1.0e-4),delzmin*(0.5-1.e-4)) !min. denominator for [a,s,t]mat 

!...  Bottom friction
      read(15,*) ntau
      if(ntau.lt.0.or.ntau.gt.2) then
      	write(11,*)'Unknown ntau', ntau
      	stop
      endif
      if(ntau==0) then
        read(15,*) Cd0
      else if(ntau==1) then
      	open(21,file='rough.bp',status='old')
      	read(21,*)
      	read(21,*) npbp
      	if(npbp.ne.np) then
          write(11,*)'rough.bp is hgrid.gr3 based'
          stop
      	endif
        do i=1,np
          read(21,*)j,xtmp,ytmp,swild(i)
      	enddo
      	do i=1,ns
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          roughness(i)=(swild(n1)+swild(n2))/2 !in meters
	  if(roughness(i)<0) then
	    write(*,*)'Negative bottom roughness',roughness(i),i
	    stop
	  endif
      	enddo
      	close(21)
      else !ntau=2
!	drag.bp is a build pt file because gredit cannot keep enough precision
!	for gredit depth
      	open(21,file='drag.bp',status='old')
      	read(21,*)
      	read(21,*) npbp
      	if(npbp.ne.np) then
          write(11,*)'drag.bp is hgrid.gr3 based'
          stop
      	endif
        do i=1,np
          read(21,*)j,xtmp,ytmp,swild(i)
      	enddo
      	do i=1,ns
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          bdragc(i)=(swild(n1)+swild(n2))/2
      	enddo
      	close(21)
      endif

!     Coriolis
      read(15,*) ncor
      if(abs(ncor)>1) then
        write(11,*)'Unknown ncor',ncor
        stop
      endif
      if(ncor==-1) then
        read(15,*) tmp
 	coricoef=2*omega*dsin(tmp/180*pi)
      else if(ncor==0) then
	read(15,*) coricoef
      else !ncor=1
      	write(*,*)'Check slam0 and sfea0 as variable Coriolis is used'
      	write(16,*)'Check slam0 and sfea0 as variable Coriolis is used'
      	open(32,file='hgrid.ll',status='old')
      	read(32,*)
      	read(32,*) !ne,np
      	do i=1,np
          read(32,*)j,xlon(i),ylat(i)
          xlon(i)=xlon(i)*deg2rad
          ylat(i)=ylat(i)*deg2rad
      	enddo !i
      	close(32)
      endif

!     Wind
      read(15,*) nws,wtiminc
      if(nws.lt.0.or.nws.gt.2) then
        write(11,*)'Unknown nws',nws
        stop
      endif
      if(nws.gt.0.and.dt.gt.wtiminc) then
      	write(11,*)'wtiminc < dt'
      	stop
      endif

      if(nws.gt.0) read(15,*) nrampwind,drampwind

!	Heat and salt conservation flags
      read(15,*) ihconsv,isconsv
      if(ihconsv.lt.0.or.ihconsv.gt.1.or.isconsv.lt.0.or.isconsv.gt.1) then
        write(11,*)'Unknown ihconsv or isconsv',ihconsv,isconsv
      	stop
      endif
      if(ihconsv+isconsv.ne.0.and.nws.ne.2) then
        write(11,*)'Heat/salt budge model must have nws=2'
      	stop
      endif
      if(ihconsv.ne.0) then
        write(16,*)'Warning: you have chosen a heat conservation model'
      	write(16,*)'which assumes start time at 0:00 PST!'
      endif

!...  Turbulence closure options
      read(15,*) itur
      if(itur.lt.-2.or.itur.gt.3) then
      	write(11,*)'Unknown turbulence closure model',itur
      	stop
      endif
      if(ibc.eq.1.and.ibtp.eq.0.and.itur.gt.0) then
        write(11,*)'Barotropic model cannot have turbulence closure'
        stop
      endif

      if(itur.eq.0) then
        read(15,*) vdiffcon,tdiffcon
        do j=1,nvrt
          do i=1,ns
            vdiff(i,j)=vdiffcon
            tdiff(i,j)=tdiffcon
          enddo
        enddo !j=1,nvrt
      else if(itur.eq.-1) then !VVD
        open(43,file='vvd.dat',status='old')
      	read(43,*) !nvrt
      	do j=1,nvrt
          read(43,*)k,vdiffcon,tdiffcon
          do i=1,mns
            vdiff(i,j)=vdiffcon
            tdiff(i,j)=tdiffcon
          enddo
      	enddo !j=1,nvrt
        close(43)
      else if(itur.eq.-2) then !HVD
        open(43,file='hvd.mom',status='old')
        open(44,file='hvd.tran',status='old')
      	read(43,*)
      	read(43,*) !nsbp
      	read(44,*)
      	read(44,*) !nsbp
      	do i=1,ns
          read(43,*)k,xtmp,ytmp,vdiffcon
          read(44,*)k,xtmp,ytmp,tdiffcon
          do j=1,nvrt
            vdiff(i,j)=vdiffcon
            tdiff(i,j)=tdiffcon
          enddo
        enddo !i=1,ns
      	close(43)
      	close(44)
      else if(itur.eq.2) then !read in P&P coefficients
        read(15,*) tdmin_pp,h1_pp,vdmax_pp1,vdmin_pp1,h2_pp,vdmax_pp2,vdmin_pp2
	if(h1_pp.ge.h2_pp) then
	  write(11,*)'h1_pp >= h2_pp in P&P'
	  stop
	endif
	if(vdmax_pp1.lt.vdmin_pp1.or.vdmax_pp2.lt.vdmin_pp2) then
	  write(11,*)'Wrong limits in P&P:',vdmax_pp1,vdmin_pp1,vdmax_pp2,vdmin_pp2
	  stop
	endif
      else if(itur.eq.3) then !read in const. (cf. Umlauf and Burchard 2003)
	read(15,*) mid,stab
        read(15,*) bgdiff,hestu_my,diffmax_est,xlmin_est,hcont_my,diffmax_sea,xlmin_sea
      	if(hestu_my>hcont_my) then
          write(11,*)'hestu_my > hcont_my'
          stop
        endif
        if(diffmax_est<bgdiff.or.diffmax_sea<bgdiff) then
          write(11,*)'Bg. diffusivity > diffmax'
          stop
        endif

!	Constants used in GLS; cpsi3 later
	cmiu0=dsqrt(0.3d0)
	a2_cm03=2/cmiu0**3
        eps_min=1.e-12

	select case(mid)
	  case('MY') 
	    rpub=0; rmub=1; rnub=1; cpsi1=0.9; cpsi2=0.5
	    q2min=5.e-6; psimin=1.e-8
	    if(stab.ne.'GA') then
	      write(11,*)'MY must use Galperins ASM:',stab
	      stop
	    endif
	  case('KL')
	    rpub=0; rmub=1; rnub=1; schk=2.44; schpsi=2.44; cpsi1=0.9; cpsi2=0.5
	    q2min=5.e-6; psimin=1.e-8
	  case('KE')
	    rpub=3; rmub=1.5; rnub=-1; schk=1; schpsi=1.3; cpsi1=1.44; cpsi2=1.92
	    q2min=7.6e-6; psimin=1.e-12
	  case('KW')
	    rpub=-1; rmub=0.5; rnub=-1; schk=2; schpsi=2; cpsi1=0.555; cpsi2=0.833
	    q2min=7.6e-6; psimin=1.e-12
	  case('UB')
	    rpub=2; rmub=1; rnub=-0.67; schk=0.8; schpsi=1.07; cpsi1=1; cpsi2=1.22
	    q2min=7.6e-6; psimin=1.e-12
	  case default
	    write(11,*)'Unknow closure:',mid
	    stop
	end select
	if(rnub==0) then
	  write(11,*)'Wrong input for rnub:',rnub
	  stop
	endif

	if(stab.ne.'GA'.and.stab.ne.'KC') then
	  write(11,*)'Unknown ASM:',stab
	  stop
	endif

!	Consts. used in Canuto's ASM (Model A)
   	ubl1=0.1070
   	ubl2=0.0032
   	ubl3=0.0864
   	ubl4=0.12
   	ubl5=11.9
   	ubl6=0.4
   	ubl7=0
   	ubl8=0.48
   	ubs0=1.5*ubl1*ubl5**2
   	ubs1=-ubl4*(ubl6+ubl7)+2*ubl4*ubl5*(ubl1-ubl2/3-ubl3)+1.5*ubl1*ubl5*ubl8
   	ubs2=-0.375*ubl1*(ubl6**2-ubl7**2)
   	ubs4=2*ubl5
   	ubs5=2*ubl4
   	ubs6=2*ubl5/3*(3*ubl3**2-ubl2**2)-0.5*ubl5*ubl1*(3*ubl3-ubl2)+0.75*ubl1*(ubl6-ubl7)
   	ubd0=3*ubl5**2
   	ubd1=ubl5*(7*ubl4+3*ubl8)
   	ubd2=ubl5**2*(3*ubl3**2-ubl2**2)-0.75*(ubl6**2-ubl7**2)
   	ubd3=ubl4*(4*ubl4+3*ubl8)
   	ubd4=ubl4*(ubl2*ubl6-3*ubl3*ubl7-ubl5*(ubl2**2-ubl3**2))+ubl5*ubl8*(3*ubl3**2-ubl2**2)
   	ubd5=0.25*(ubl2**2-3*ubl3**2)*(ubl6**2-ubl7**2)
!	print*, 'ubd2=',ubd2,',ubd4=',ubd4,',ubd2/ubd4=',ubd2/ubd4

        do i=1,ns
	  xlmin2(i)=2*q2min*0.1*dmax1(h0,dps(i)) !floor for non-surface layers
          if(dps(i)<=hestu_my) then
            xlmin1(i)=xlmin_est
	    diffmax(i)=diffmax_est
          else if(dps(i)<=hcont_my) then
            xlmin1(i)=xlmin_est+(xlmin_sea-xlmin_est)*(dps(i)-hestu_my)/(hcont_my-hestu_my)
	    diffmax(i)=diffmax_est+(diffmax_sea-diffmax_est)*(dps(i)-hestu_my)/(hcont_my-hestu_my)
          else !dps > hcont_my
            xlmin1(i)=xlmin_sea
	    diffmax(i)=diffmax_sea
          endif
          xlmin1(i)=dmax1(xlmin1(i),xlmin2(i))
        enddo !i

      endif !itur

!...  HORCON parameter
      read(15,*) ihorcon,smagcoef

!      if(ihorcon.eq.0) then
!	read(15,*) hor 
!	do i=1,mns
!	  smagcoef(i)=hor
!	enddo
!      else !ihorcon=1
!	open(99,file='fort.99',status='old')
!	read(99,*)
!	read(99,*) nsbp
!	do i=1,nsbp
!	  read(99,*)j,xtmp,ytmp,hor
!	  smagcoef(i)=hor
!	enddo
!	close(99)
!      endif

      read(15,*) ihorcon_st
      if(ihorcon_st.eq.0) then
        read(15,*) hor_t0,hor_s0 
        do i=1,ne
          hor_t(i)=hor_t0
          hor_s(i)=hor_s0
        enddo
      else !ihorcon_st=1
      	open(97,file='fort.99_T',status='old')
      	open(98,file='fort.99_S',status='old')
      	read(97,*)
      	read(97,*) !ne
      	read(98,*)
      	read(98,*) !ne
      	do i=1,ne
          read(97,*)j,xtmp,ytmp,hor_t0
          read(98,*)j,xtmp,ytmp,hor_s0
          hor_t(i)=hor_t0
          hor_s(i)=hor_s0
        enddo
        close(97)
        close(98)
      endif

      read(15,*)ictemp,icsalt
      if(ictemp.ne.1.and.ictemp.ne.2.or.icsalt.ne.1.and.icsalt.ne.2) then
      	write(11,*)'Unknown i.c. flag',ictemp,icsalt
      	stop
      endif

!...  Sponge layer
      read(15,*) isponge
      if(isponge.ne.0) then
        open(96,file='sponge.gr3',status='old')
      	read(96,*)
      	read(96,*) !ne,np
      	do i=1,np
          read(96,*)j,xtmp,ytmp,relax(i)
      	enddo !i
      	close(96)
      endif

!...  Earth tidal potential
      read(15,*) ntip,tip_dp !cut-off depth for applying tidal potential
      if(ntip.gt.mnbfr) then
        write(11,*)'ntip > mnbfr',ntip,mnbfr
        stop
      endif
      
      if(ntip.gt.0) then
	open(32,file='hgrid.ll',status='old')
        read(32,*)
        read(32,*) !ne,np
        do i=1,np
          read(32,*)j,xlon(i),ylat(i)
          xlon(i)=xlon(i)*deg2rad
          ylat(i)=ylat(i)*deg2rad
!...      Pre-compute species function to save time
	  fun_lat(i,0)=3*dsin(ylat(i))**2-1
	  fun_lat(i,1)=dsin(2*ylat(i))
	  fun_lat(i,2)=dcos(ylat(i))**2
        enddo !i
        close(32)

	do i=1,ne
	  xlon_e(i)=0
	  ylat_e(i)=0
	  do j=1,i34(i)
	    xlon_e(i)=xlon_e(i)+xlon(nm(i,j))/i34(i)
	    ylat_e(i)=ylat_e(i)+ylat(nm(i,j))/i34(i)
	  enddo !j
	  fun_lat_e(i,0)=3*dsin(ylat_e(i))**2-1
	  fun_lat_e(i,1)=dsin(2*ylat_e(i))
	  fun_lat_e(i,2)=dcos(ylat_e(i))**2
	enddo !i
      endif !ntip.gt.0
	
      do i=1,ntip
        read(15,*) !tag
	read(15,*) jspc(i),tamp(i),tfreq(i),tnf(i),tear(i)
	if(jspc(i).lt.0.or.jspc(i).gt.2) then
	  write(11,*)'Illegal tidal species #',jspc(i)
	  stop
	endif
        tear(i)=tear(i)*deg2rad
      enddo !i

!...  Boundary forcing freqs.
      read(15,*) nbfr
      if(nbfr.gt.mnbfr) then
        write(11,*)'nbfr > mnbfr',nbfr,mnbfr
        stop
      endif

      do i=1,nbfr
        read(15,'(a5)') bountag(i)
        read(15,*) amig(i),ff(i),face(i) !freq., nodal factor and earth equil.
        face(i)=face(i)*deg2rad
      enddo

      read(15,*) nope1
      if(nope1.ne.nope) then
        write(11,*)'Inconsistent # of open bnds',nope1,nope
      	stop
      endif

      nettype=0 !total # of type I eta bnds
      nfltype=0
      ntetype=0
      nsatype=0
      do k=1,nope
        read(15,*) ntmp,iettype(k),ifltype(k),itetype(k),isatype(k)
      	if(ntmp.ne.noe(k)) then
          write(11,*)'Inconsistent # of elements at open boundary',k
          write(11,*)ntmp,noe(k)
          stop
        endif

        if(iettype(k).eq.1) then
          nettype=nettype+1
!	  Mock reading
          open(50,file='elev.th',status='old')
          do j=1,nt
            read(50,*) ttt,et
          enddo !j
          rewind(50)
        else if(iettype(k).eq.2) then
          read(15,*) eth(k)
        else if(iettype(k).eq.3) then
          do i=1,nbfr
            read(15,'(a10)') fq_nm(i) !freq. name
            do j=1,noe(k)
              read(15,*) emo(k,i,j),efa(k,i,j) !amp. and phase
              efa(k,i,j)=efa(k,i,j)*deg2rad
            enddo
          enddo
        else if(iettype(k).ne.-1.and.iettype(k).ne.0) then
          write(11,*)'INVALID VALUE FOR IETTYPE'
          stop
        endif

        if(iabs(ifltype(k))==1) then
          nfltype=nfltype+1
          open(51,file='flux.th',status='old')
          do j=1,nt
            read(51,*) ttt,qq
          enddo
          rewind(51)
      	else if(iabs(ifltype(k))==2) then
          read(15,*) qth(k)
        else if(ifltype(k)/=0) then
          write(11,*) 'INVALID VALUE FOR IFLTYPE'
          stop
        endif

	if(ifltype(k)<0.and.iettype(k)/=0) then
	  write(11,*)'Flather b.c. requires iettype=0'
	  stop
	endif

        if(itetype(k).eq.1) then
          ntetype=ntetype+1
          open(52,file='temp.th',status='old')
          do j=1,nt
            read(52,*) ttt,temp
          enddo
          rewind(52)
        else if(itetype(k).eq.2) then
          read(15,*) tth(k)
      	else if(itetype(k).eq.-1) then
          read(15,*) tth(k)
        else if(itetype(k).lt.-1.or.itetype(k).gt.3) then
          write(11,*) 'INVALID VALUE FOR ITETYPE'
          stop
        endif

        if(isatype(k).eq.1) then
          nsatype=nsatype+1
          open(53,file='salt.th',status='old')
          do j=1,nt
            read(53,*) ttt,sal
          enddo
          rewind(53)
        else if(isatype(k).eq.2) then
          read(15,*) sth(k)
        else if(isatype(k).eq.-1) then
          read(15,*) sth(k)
        else if(isatype(k).lt.-1.or.isatype(k).gt.3) then
          write(11,*) 'INVALID VALUE FOR ISATYPE'
          stop
        endif
      enddo !k=1,nope

!     Flow sides 
      do i=1,ns
      	isflowside(i)=0
      	isflowside2(i)=0
      enddo

      do i=1,nope
        if(ifltype(i)>=1) then 
          do j=1,noe(i)
            iel=ioe(i,j)
            do k=1,i34(iel)
              isd=js(iel,k)
              isflowside(isd)=i
            enddo !k
      	  enddo !j
        endif

	if(ifltype(i)<=-1) then
	  do j=1,nosd(i)
	    isd=iosd(i,j)
	    isflowside2(isd)=i
	    if(dps(isd)<=0) then
	      write(11,*)'Depth < 0 at radiation bnd side:',i,j,isd
	      stop
	    endif
	  enddo
	endif
      enddo

!     Global output parameters
      noutput=20
      if(noutput.gt.mnout) then
        write(11,*)'Increase mnout in the header to',noutput
        stop
      endif
      outfile(1)='elev.61'
      outfile(2)='pres.61'
      outfile(3)='airt.61'
      outfile(4)='shum.61'
      outfile(5)='srad.61'
      outfile(6)='flsu.61'
      outfile(7)='fllu.61'
      outfile(8)='radu.61'
      outfile(9)='radd.61'
      outfile(10)='flux.61'
      outfile(11)='wind.62'
      outfile(12)='wist.62'
      outfile(13)='hvel.64'
      outfile(14)='vert.63'
      outfile(15)='temp.63'
      outfile(16)='salt.63'
      outfile(17)='conc.63'
      outfile(18)='tdff.63'
      outfile(19)='kine.63'
      outfile(20)='mixl.63'
      variable_nm(1)='surface elevation'
      variable_nm(2)='atmopheric pressure'
      variable_nm(3)='air temperature'
      variable_nm(4)='specific humidity'
      variable_nm(5)='solar radiation'
      variable_nm(6)='fluxsu'
      variable_nm(7)='fluxlu'
      variable_nm(8)='hradu'
      variable_nm(9)='hradd'
      variable_nm(10)='total flux'
      variable_nm(11)='wind speed'
      variable_nm(12)='wind stress (m^2/s^2)'
      variable_nm(13)='horizontal velocity'
      variable_nm(14)='vertical velocity'
      variable_nm(15)='temperature in C'
      variable_nm(16)='salinity in psu'
      variable_nm(17)='density in kg/m^3'
      variable_nm(18)='diffusivity for transport'
      variable_nm(19)='turbulent kinetic energy'
      variable_nm(20)='turbulent mixing length'
      do i=1,10
        variable_dim(i)='2D scalar'
      enddo
      variable_dim(11)='2D vector'
      variable_dim(12)='2D vector'
      variable_dim(13)='3D vector'
      do i=14,noutput
        variable_dim(i)='3D scalar'
      enddo
      do i=1,noutput
        if(i.le.12) then
          vpos(i)=0
      	else if(i.ge.15.and.i.le.17) then
          vpos(i)=0.5
      	else 
          vpos(i)=1
      	endif
      enddo

      read(15,*) nspool,ihfskip !output and file spools
      if(nspool.eq.0.or.ihfskip.eq.0) then
        write(11,*)'Zero nspool'
        stop
      endif
      if(mod(ihfskip,nspool).ne.0) then
	write(11,*)'ihfskip/nspool != integer'
	stop
      endif
      nrec=min(nt,ihfskip)/nspool

      do i=1,noutput
        read(15,*) iof(i)
        if(iof(i).ne.0.and.iof(i).ne.1) then
          write(11,*)'Unknown output option',i,iof(i)
          stop
        endif
      enddo !i=1,noutput

!...  serg
!...  read from the parameter file output options for da varaibles
      da_noutput = 4
      da_outfile(1)='eta2.61'
      da_outfile(2)='tsd2.63'
      da_outfile(3)='ssd2.63'
      da_outfile(4)='usg_.64'
      da_variable_nm(1)='elevation at elem'
      da_variable_nm(2)='temperature at side centers'
      da_variable_nm(3)='salinity at side centers'
      da_variable_nm(4)='velocity in xy at side centers'
      da_variable_dim(1)='2D scalar'
      da_variable_dim(2)='3D scalar'
      da_variable_dim(3)='3D scalar'
      da_variable_dim(4)='3D vector'
      da_vpos(1)=0
      da_vpos(2)=0.5
      da_vpos(3)=0.5
      da_vpos(4)=0.5

      do i=1,da_noutput
        read(15,*) da_iof(i)
        if(iof(i).ne.0.and.iof(i).ne.1) then
          write(11,*)'Unknown output option',i,da_iof(i)
          stop
        endif
      enddo !i=1,da_noutput
!...  serg end

!...  Test output parameters
      read(15,*) noutgm 
      if(noutgm.ne.1.and.noutgm.ne.0) then
        write(11,*)'Unknown noutgm',noutgm
        stop
      endif
      
!...  input information about hot start output
!...
      read(15,*) nhstar
      if(nhstar.ne.0.and.nhstar.ne.1) then
        write(11,*)'Unknown nhstar',nhstar
        stop
      endif

!... Itpack solver info
      read(15,*) isolver,itmax1,iremove,zeta,tol
      if(itmax1.gt.itmax) then
        write(11,*)'Increase itmax in header file'
      	stop
      endif
      if(isolver.lt.1.or.isolver.gt.4) then
        write(11,*)'Unknown solver',isolver
        stop
      endif

!...  Compute flux flag
      read(15,*) iflux,ihcheck

!.... W switch; inactive
      read(15,*) !iwmode

!.... Mode splitting factor
      read(15,*) nsplit

!...  Eqstate
      read(15,*) ieqstate
      if(ieqstate<0.or.ieqstate>1) then
	write(11,*)'Unknown ieqstate:',ieqstate
	stop
      endif

!...  Check last parameter read in from fort.15
      write(*,*)'Last parameter in param.in is ieqstate=',ieqstate
      close(15)
!     End reading fort.15

      if(nscreen.eq.1) write(*,*)'done reading fort.14, 15, and 19...'
      write(16,*)'done reading fort.14, 15, and 19...'


!								   *
!*******************************************************************
!								   *
!	Initialization for cold and hot start			   *
!								   *
!*******************************************************************
!								   *

!...  Coriolis parameter
!...
      if(ncor<=0) then
      	do i=1,ns
          cori(i)=coricoef
      	enddo
      else !ncor=1
        open(31,file='coriolis.out')
      	fc=2*omega*dsin(sfea0)
      	beta=2*omega*dcos(sfea0)
      	do i=1,ns
          id1=isidenode(i,1)
          id2=isidenode(i,2)
          sphi=(ylat(id1)+ylat(id2))/2
          cori(i)=fc+beta*(sphi-sfea0)
	  if(iwrite.eq.0) then
            write(31,*)i,xcj(i),ycj(i),cori(i)
	  else !evm
	    write(31,"(a,i6,a,f16.9,a,f16.9,a,es22.14e3,a)",advance="no") &
     &           " ",i," ",xcj(i)," ",ycj(i)," ",cori(i),"\n"
	  endif
      	enddo !i=1,ns
      	close(31)
      endif


!								   *
!*******************************************************************
!								   *
!	Initialization for cold start alone			   *
!								   *
!*******************************************************************
!								   *

      if(ihot.eq.0) then
!------------------------------------------------------------------
!...  read the initial salinity and temperature values from 
!...  salinity.bp and temperature.bp files. Initial S,T fields may vary
!...  either horizontally (and vertically homogeneous) or vertically 
!...  (horizontally homogeneous). For more general 3D case, use hot start.
!...
      if(ibc.eq.1.and.ibtp.eq.0) then
!	Reset ictemp and icsalt
      	ictemp=1
      	icsalt=1
        if(10.lt.tempmin.or.10.gt.tempmax.or.33.lt.saltmin.or.33.gt.saltmax) then
          write(11,*)'Pls reset ST range to include S=33 and T=10C'
          stop
        endif
      	do i=1,np
	  do k=1,nvrt
            tem0(i,k)=10
            sal0(i,k)=33
	  enddo !k
      	enddo !i
      else !read in S,T
        open(24,file='temp.ic',status='old')
        open(25,file='salt.ic',status='old')
        if(ictemp.eq.1) then
          read(24,*) 
          read(24,*) !np
          do i=1,np
            read(24,*) num,xtmp,ytmp,te
            if(te.lt.tempmin.or.te.gt.tempmax) then
              write(11,*)'Initial invalid T at',i,te
              stop
            endif
	    do k=1,nvrt
              tem0(i,k)=te
	    enddo !k
          enddo !i
        else !ictemp=2
          read(24,*) !nvrt
          do k=1,nvrt
            read(24,*)num,te
            if(te.lt.tempmin.or.te.gt.tempmax) then
              write(11,*)'Initial invalid T at',k,te
              stop
            endif
	    do i=1,np
              tem0(i,k)=te
	    enddo !i
          enddo !k
        endif

        if(icsalt.eq.1) then
          read(25,*) 
          read(25,*) !np
          do i=1,np
            read(25,*) num,xtmp,ytmp,sa
            if(sa.lt.saltmin.or.sa.gt.saltmax) then
              write(11,*)'Initial invalid S at',i,sa
              stop
            endif
	    do k=1,nvrt
              sal0(i,k)=sa
	    enddo !k
          enddo
        else !icsalt=2
          read(25,*) !nvrt
          do k=1,nvrt
            read(25,*)num,sa
            if(sa.lt.saltmin.or.sa.gt.saltmax) then
              write(11,*)'Initial invalid S at',k,sa
              stop
            endif
	    do i=1,np
              sal0(i,k)=sa
	    enddo !i
          enddo !i
        endif
        close(24)
        close(25)
      endif !ibc.eq.1.and.ibtp.eq.0

!...  initialize S,T
      do i=1,np
        do j=1,nvrt
          tnd(i,j)=tem0(i,j)
          snd(i,j)=sal0(i,j)
      	enddo
      enddo

      do i=1,ns
	n1=isidenode(i,1)
	n2=isidenode(i,2)
        do j=1,nvrt
          tsd(i,j)=(tem0(n1,j)+tem0(n2,j))/2
          ssd(i,j)=(sal0(n1,j)+sal0(n2,j))/2
      	enddo
      enddo

!...  initialize elevations and vel. 
!...
      do i=1,ne
        eta1(i)=0
        eta2(i)=0
        do j=0,nvrt
          we(i,j)=0
        enddo
      enddo

      do i=1,np
        peta(i)=0
      	ibad(i)=0
        do k=0,nvrt
          uu1(i,k)=0.!only for internal mode
          vv1(i,k)=0
          ww1(i,k)=0
      	enddo
      enddo 

      do i=1,ns
        do j=0,nvrt
          vn2(i,j)=0 !0.1*snx(i)
          vt2(i,j)=0 !-0.1*sny(i)
        enddo
      enddo

!... initialize q2 and xl in MY-G scheme
!...
      if(itur.eq.3) then
        do i=1,ns
      	  do j=0,nvrt
            xl(i,j)=xlmin2(i)
            q2(i,j)=q2min
          enddo !j
        enddo !i=1,ns	
      endif !itur=3

!...  initialize wind for nws=1,2 (first two lines)
!...
      if(nws.eq.1) then
        open(22,file='wind.th',status='old')
        read(22,*) wx1,wy1
        read(22,*) wx2,wy2
        do i=1,np
          windx1(i)=wx1
          windy1(i)=wy1
          windx2(i)=wx2
          windy2(i)=wy2
        enddo
        wtime1=0
      	wtime2=wtiminc 
      endif

!	CORIE mode
      if(nws.eq.2) then
      	wtime1=0
      	wtime2=wtiminc 
      	call get_wind(wtime1,windx1,windy1,pr1,airt1,shum1)
      	call get_wind(wtime2,windx2,windy2,pr2,airt2,shum2)
      endif !nws=2
!------------------------------------------------------------------
      endif !ihot=0

      if(nscreen.eq.1) write(*,*)'done initializing cold start'
      write(16,*)'done initializing cold start'

!                                                                             
!******************************************************************************
!                                                                             *
!		hot start setup of the program				      *
!                                                                             *
!******************************************************************************
!
!	Record length for hot start files (double precision for all reals)
      ihot_len=nbyte*(3+4*ne+2*ne*(nvrt+1)+4*ns*(nvrt+1)+4*ns*nvrt+ &
     &	       3*np+7*np*(nvrt+1)+8*np*nvrt+1)+12
      if(itur.eq.3) ihot_len=ihot_len+nbyte*4*ns*(nvrt+1)

      if(ihot.ne.0) then
        if (da_ihot_type.eq.1) then !serg hotstart in elcirc format
          print*, 'starting from elcirc-type hotstart'
          open(36,file='hotstart.in',access='direct',recl=ihot_len)
          if(itur==3) then
            read(36,rec=1)time,iths,(eta1(i),eta2(i), &
            &(we(i,j),j=0,nvrt),i=1,ne), &
            &((vn2(i,j),vt2(i,j),j=0,nvrt),(tsd(i,j),ssd(i,j),j=1,nvrt),i=1,ns) &
            &,(peta(i),ibad(i),(uu1(i,j),vv1(i,j),ww1(i,j),junk,j=0,nvrt), &
            &(tnd(i,j),snd(i,j),tem0(i,j),sal0(i,j),j=1,nvrt),i=1,np), &
            &((q2(i,j),xl(i,j),j=0,nvrt),i=1,ns), &
            &ifile,ifile_char

            do i=1,ns
              do j=0,nvrt
                q2(i,j)=dmax1(q2min,q2(i,j))
                xl(i,j)=dmax1(xlmin2(i),xl(i,j))
              enddo
            enddo
          else !itur.ne.3
            read(36,rec=1)time,iths,(eta1(i),eta2(i), &
            &(we(i,j),j=0,nvrt),i=1,ne), &
            &((vn2(i,j),vt2(i,j),j=0,nvrt),(tsd(i,j),ssd(i,j),j=1,nvrt),i=1,ns) &
            &,(peta(i),ibad(i),(uu1(i,j),vv1(i,j),ww1(i,j),junk,j=0,nvrt), &
            &(tnd(i,j),snd(i,j),tem0(i,j),sal0(i,j),j=1,nvrt),i=1,np), &
            &ifile,ifile_char
          endif !itur
          close(36)
        else if(da_ihot_type.eq.2) then !serg hotstart in da format
          print*, 'starting from da-type hotstart'
          ihot_len_da=nbyte*(3+4*ne+2*ne*(nvrt+1)+4*ns*(nvrt+1)+4*ns*nvrt+ &
          &         8*np*nvrt+1)+12
          if(itur.eq.3) ihot_len_da=ihot_len_da+nbyte*4*ns*(nvrt+1)

          open(36,file='hotstart.in',access='direct',recl=ihot_len_da)
          if(itur==3) then
            read(36,rec=1)time,iths,(eta1(i),eta2(i), &
            &(we(i,j),j=0,nvrt),i=1,ne), &
            &((da_ugsc(i,j),da_vgsc(i,j),j=0,nvrt),(tsd(i,j),ssd(i,j),j=1,nvrt),i=1,ns) &
            &,((tnd(i,j),snd(i,j),tem0(i,j),sal0(i,j),j=1,nvrt),i=1,np), &
            &((q2(i,j),xl(i,j),j=0,nvrt),i=1,ns), &
            &ifile,ifile_char

            do i=1,ns
              do j=0,nvrt
                q2(i,j)=dmax1(q2min,q2(i,j))
                xl(i,j)=dmax1(xlmin2(i),xl(i,j))
              enddo
            enddo
          else !itur.ne.3
            read(36,rec=1)time,iths,(eta1(i),eta2(i), &
            &(we(i,j),j=0,nvrt),i=1,ne), &
            &((da_ugsc(i,j),da_vgsc(i,j),j=0,nvrt),(tsd(i,j),ssd(i,j),j=1,nvrt),i=1,ns) &
            &,((tnd(i,j),snd(i,j),tem0(i,j),sal0(i,j),j=1,nvrt),i=1,np), &
            &ifile,ifile_char
          endif !itur
          close(36)

          do i=1,ne
            eta1(i)=eta2(i)
          enddo

!          do i=1,np
!            peta(i)=0
!            ibad(i)=0
!            do k=0,nvrt
!              uu1(i,k)=0.!only for internal mode
!              vv1(i,k)=0
!              ww1(i,k)=0
!            enddo
!          enddo

          do i=1,ns
            do j=0,nvrt
              vn2(i,j)=  da_ugsc(i,j)*snx(i)+da_vgsc(i,j)*sny(i)
              vt2(i,j)= -da_ugsc(i,j)*sny(i)+da_vgsc(i,j)*snx(i)
            enddo
          enddo
        endif !serg endif da_ihot_type 
!...    change time and iteration for forecast mode
!...    Causion: this affects all t.h. files (fort.5[0-3]) and wind files
        if(ihot==1) then
          time=0
          iths=0
        endif

        write(*,*)'hot start at time=',time,iths
        write(16,*)'hot start at time=',time,iths

!...  find position in the wind input file for nws=1,2, and read in wind[x,y][1,2]
!...
        if(nws.eq.1) then
          open(22,file='wind.th',status='old')
          rewind(22)
          ninv=time/wtiminc
          wtime1=ninv*wtiminc 
          wtime2=(ninv+1)*wtiminc 
          do it=0,ninv
            read(22,*)wx1,wy1
          enddo
          read(22,*)wx2,wy2
          do i=1,np
            windx1(i)=wx1
            windy1(i)=wy1
            windx2(i)=wx2
            windy2(i)=wy2
          enddo
        endif

        if(nws.eq.2) then
          ninv=time/wtiminc
          wtime1=ninv*wtiminc 
          wtime2=(ninv+1)*wtiminc 
          call get_wind(wtime1,windx1,windy1,pr1,airt1,shum1)
          call get_wind(wtime2,windx2,windy2,pr2,airt2,shum2)
        endif !nws=2

!...  Find positions in t.h. files fort.5[0-3] 
        if(nettype>0) then
          do it=1,iths
            read(50,*) ttt,et
          enddo !it
        endif

        if(nfltype>0) then
          do it=1,iths
            read(51,*) ttt,qq
          enddo !it
        endif

      	if(ntetype>0) then
          do it=1,iths
            read(52,*) ttt,te
          enddo !it
      	endif

      	if(nsatype>0) then
          do it=1,iths
            read(53,*) ttt,sal
          enddo !it
      	endif

!...  end hot start section
!...
      endif !ihot.ne.0

      if(nscreen.eq.1) write(*,*)'done initializing variables...'
      write(16,*)'done initializing variables...'

!                                                                             *
!******************************************************************************
!                                                                             *
!			open output files				      *
!                                                                             *
!******************************************************************************
!                                                                             *

!...  write global output headers
!...

      if(ihot<=1) then
   	ifile=1 !output file #
!       Convert it to a string
      	write(ifile_char,'(i12)') ifile
      endif
      call header
!... serg
      call da_header
!... serg

      if(nscreen.eq.1) write(*,*)'done initializing outputs'
      write(16,*)'done initializing outputs'

!                                                                             *
!******************************************************************************
!                                                                             *
!	    	          Time stepping 				      *
!                                                                             *
!******************************************************************************
!                                                                             *

      if(ihot.eq.0) iths=0
      if(nscreen.eq.1) write(*,*)'time stepping begins...',iths+1,nt
      write(16,*)'time stepping begins...',iths+1,nt

!...  Compute initial levels
      call levels(iths,iths,itur)
      if(nscreen.eq.1) write(*,*)'done computing initial levels...'
      write(16,*)'done computing initial levels...'

!... serg compute peta and ibad
      if (da_ihot_type.eq.2) then
!...  Compute nodal elevations for tang. vel.
        do i=1,np
          ibad(i)=0 !reliable elevation

          sum1=0 !sum of elemental elevations
          sum=0 !sum of areas
          do j=1,nne(i)
            ie=ine(i,j)
            if(kfe(ie).eq.0) then
              ibad(i)=1
            else
              sum1=sum1+areas(ie)*eta2(ie)
              sum=sum+areas(ie)
            endif
          enddo !j=1,nne(i)

          if(ibad(i).eq.0) peta(i)=sum1/sum
        enddo !i=1,np

!     Estimate bad elevations for rewetting
        do i=1,np
          if(ibad(i).eq.1) then
            sum=0
            icount=0
            do j=1,nnp(i)
              nd=inp(i,j)
              if(ibad(nd).eq.0) then
                icount=icount+1
                sum=sum+peta(nd)
              endif
            enddo !j
            if(icount.ne.0) peta(i)=sum/icount
          endif !ibad(i).eq.1
        enddo !i=1,np
      endif
!... serg

!...  Compute total # of active faces
      nfaces=0
      do i=1,ns
        if(kfs(i)/=0) nfaces=nfaces+kfs(i)-kbs(i)+1
      enddo !i=1,ns
      if(nscreen.eq.1) write(*,*)'Total # of faces=',nfaces
      write(16,*)'Total # of faces=',nfaces

!...  Compute nodal vel. 
      call nodalvel
      if(nscreen.eq.1) write(*,*)'done computing initial nodal vel...'
      write(16,*)'done computing initial nodal vel...'
!	serg assign values of uu ...
      if (da_ihot_type.eq.2) then
          uu1=uu2
          vv1=vv2
          ww1=ww2
      endif
!     serg
!...  Compute initial density
      call eqstate
      if(nscreen.eq.1) write(*,*)'done computing initial density...'
      write(16,*)'done computing initial density...'

!...  Store initial elevation field (hot start) for relaxation in sponge layer
!...
      do i=1,ne
        etaic(i)=eta1(i)
      enddo

!...  For UB & MYG, compute vertical diffusivities at whole levels 
      if(itur.eq.3) then
        do i=1,ns
    	  do j=kbs(i),kfs(i) !wet sides
	    call asm(g,i,j,vd,td,qd,qd2)
            vdiff(i,j)=dmin1(diffmax(i),dmax1(bgdiff,vd))
            tdiff(i,j)=dmin1(diffmax(i),dmax1(bgdiff,td))
            qdiff(i,j)=dmin1(diffmax(i),dmax1(bgdiff,qd))
            qdiff2(i,j)=dmin1(diffmax(i),dmax1(bgdiff,qd2))
          enddo !j=kbs(i),kfs(i)
        enddo !i=1,ns
      endif !itur=3

!...  Initialize heat budget model
      if(nws==2.and.ihconsv/=0) then
        call surf_fluxes(wtime1,windx1,windy1,pr1,airt1,shum1,&
     &srad,fluxsu,fluxlu,hradu,hradd,tauxz,tauyz)
      	do i=1,np
          sflux(i)=-fluxsu(i)-fluxlu(i)-(hradu(i)-hradd(i))
      	enddo
      	if(nscreen.eq.1) write(*,*)'heat budge model completes...'
      	write(16,*)'heat budge model completes...'
      endif

!...
!...  Begin time stepping
!...
      do it=iths+1,nt

      time=it*dt 

!...  define ramp function for boundary elevation forcing, wind and pressure
!...  forcing and tidal potential forcing
!...
      if(ibc.eq.0.and.nrampbc.ne.0) then
        rampbc=tanh(2*time/86400/drampbc)
      else
      	rampbc=1
      endif

      if(nws.gt.0.and.nrampwind.ne.0) then
        rampwind=tanh(2*time/86400/drampwind)
      else
        rampwind=1
      endif

      if(nramp.eq.1) then
 	ramp=tanh(2*time/86400/dramp)
      	ramptide=ramp
      else
        ramp=1
        ramptide=1
      endif

!...  Bottom friction apart from the vel. magnitude
!	Error: shall we change tk back to no-btrack?
      do i=1,ns
      	if(ntau==0) then
          Cd(i)=Cd0
      	else if(ntau==1) then
          if(dps(i).le.hestu) then
            bthick=bthick_estu !thickness of boundary layer (in meters)
          else if(dps(i).le.hcont) then
            bthick=bthick_cont
          else !dps(i) > hcont
            bthick=bthick_ocea
          endif
!          gthick=dmax1(dz(i,kbs(i))/2,bthick_estu) !thickness of half bottom layer
          gthick=dz(i,kbs(i))/2 !thickness of half bottom layer
	  if(roughness(i)<1.e-10) then
            Cd(i)=0
	  else
            Cd(i)=1/(2.5*dlog(gthick/roughness(i)))**2 
	  endif
!          Cdmin=1/(2.5*dlog(bthick/1.0d-2))**2 
!          Cd(i)=dmax1(Cd(i),Cdmin)
        else !ntau=2
          Cd(i)=bdragc(i)
        endif
        tk(i)=0 !initialize for some weird sides
      enddo !i=1,ns

      if(nscreen.eq.1) write(*,*)'done Smag. and bottom fric...'
      write(16,*)'done Smag. and bottom fric...'
      
!...  process new wind info 
!...
      if(nws.eq.1) then
        if(time.gt.wtime2) then
          wtime1=wtime2
          wtime2=wtime2+wtiminc
          read(22,*) wx2,wy2
          do i=1,np
            windx1(i)=windx2(i)
            windy1(i)=windy2(i)
            windx2(i)=wx2
            windy2(i)=wy2
          enddo
        endif

        wtratio=(time-wtime1)/wtiminc
!	Compute wind at nodes for output only
        do i=1,np
          windxp(i)=windx1(i)+wtratio*(windx2(i)-windx1(i))
          windyp(i)=windy1(i)+wtratio*(windy2(i)-windy1(i))
        enddo !i

        do i=1,ns
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          windx1_s=(windx1(n1)+windx1(n2))/2
          windx2_s=(windx2(n1)+windx2(n2))/2
          windy1_s=(windy1(n1)+windy1(n2))/2
          windy2_s=(windy2(n1)+windy2(n2))/2
          windx(i)=windx1_s+wtratio*(windx2_s-windx1_s)
          windy(i)=windy1_s+wtratio*(windy2_s-windy1_s)
        enddo
      endif !nws=1

!	CORIE mode
      if(nws.eq.2) then
        if(time>=wtime2) then
!...      Heat budget & wind stresses
          if(ihconsv.ne.0) then
            call surf_fluxes(wtime2,windx2,windy2,pr2,airt2,shum2,&
     &srad,fluxsu,fluxlu,hradu,hradd,tauxz,tauyz)
            do i=1,np
     	      sflux(i)=-fluxsu(i)-fluxlu(i)-(hradu(i)-hradd(i))
            enddo
            if(nscreen.eq.1) write(*,*)'heat budge model completes...'
            write(16,*)'heat budge model completes...'
          endif !ihconsv.ne.0

          wtime1=wtime2
          wtime2=wtime2+wtiminc
          do i=1,np
            windx1(i)=windx2(i)
            windy1(i)=windy2(i)
	    pr1(i)=pr2(i)
            airt1(i)=airt2(i)
            shum1(i)=shum2(i)
          enddo
          call get_wind(wtime2,windx2,windy2,pr2,airt2,shum2)

        endif !time.ge.wtime2

        wtratio=(time-wtime1)/wtiminc

!	Compute wind at nodes for output only
      	do i=1,np
          windxp(i)=windx1(i)+wtratio*(windx2(i)-windx1(i))
          windyp(i)=windy1(i)+wtratio*(windy2(i)-windy1(i))
        enddo !i

        do i=1,ns
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          windx1_s=(windx1(n1)+windx1(n2))/2
          windx2_s=(windx2(n1)+windx2(n2))/2
          windy1_s=(windy1(n1)+windy1(n2))/2
          windy2_s=(windy2(n1)+windy2(n2))/2
          windx(i)=windx1_s+wtratio*(windx2_s-windx1_s)
          windy(i)=windy1_s+wtratio*(windy2_s-windy1_s)
        enddo
      endif !nws=2

!...  compute wind stress components
      dragcmin=1.0d-3*(0.61+0.063*6)
      dragcmax=1.0d-3*(0.61+0.063*50)
      do i=1,ns
        if(nws.eq.0) then
	  tau_n(i)=0
	  tau_t(i)=0
	else if(nws.eq.1.or.nws.eq.2.and.ihconsv.eq.0) then
          wmag=dsqrt(windx(i)**2+windy(i)**2)
          windn=windx(i)*snx(i)+windy(i)*sny(i)
          windt=-windx(i)*sny(i)+windy(i)*snx(i)
          dragcoef=1.0d-3*(0.61+0.063*wmag)
          dragcoef=dmin1(dmax1(dragcoef,dragcmin),dragcmax)
          tau_n(i)=dragcoef*0.001293*wmag*windn*rampwind
          tau_t(i)=dragcoef*0.001293*wmag*windt*rampwind
	else !nws=2 and ihconsv !=0; tauxz and tauyz defined
	  n1=isidenode(i,1)
	  n2=isidenode(i,2)
	  if(kfp(n1).eq.0.or.kfp(n2).eq.0) then
	    tau_n(i)=0
	    tau_t(i)=0
	  else
	    tau_n(i)=(tauxz(n1)+tauxz(n2))/2*snx(i)+(tauyz(n1)+tauyz(n2))/2*sny(i)
	    tau_t(i)=-(tauxz(n1)+tauxz(n2))/2*sny(i)+(tauyz(n1)+tauyz(n2))/2*snx(i)
	    tau_n(i)=-tau_n(i)/rho0*rampwind !sign and scale difference between stresses tauxz and tau_n
	    tau_t(i)=-tau_t(i)/rho0*rampwind
	  endif
	endif !nws
      enddo !i=1,ns


      if(nscreen.eq.1) write(*,*)'done adjusting wind stress ...'
      write(16,*)'done adjusting wind stress ...'

!...  Get new t.h. values fort.5[0-3]
!...
      if(nettype.gt.0) then
        read(50,*) ttt,(ath(i),i=1,nettype)
        if(it.eq.iths+1.and.abs(ttt-time).gt.1.e-4) then
          write(11,*)'Starting time wrong for eta',it,ttt
          stop
        endif
      
        icount=0
        do k=1,nope
          if(iettype(k).eq.1) then
            icount=icount+1
            if(icount.gt.nettype) then
              write(11,*)'Wrong counting 1'
              stop
            endif
            eth(k)=ath(icount)
          endif
        enddo 
      endif !nettype.gt.0

      if(nfltype.gt.0) then
        read(51,*) ttt,(ath(i),i=1,nfltype)
        if(it.eq.iths+1.and.abs(ttt-time).gt.1.e-4) then
          write(11,*)'Starting time wrong for flux',it,ttt
          stop
        endif

        icount=0
      	do k=1,nope
          if(iabs(ifltype(k))==1) then
            icount=icount+1
            if(icount.gt.nfltype) then
              write(11,*)'Wrong counting 2'
              stop
            endif
            qth(k)=ath(icount)
          endif
        enddo !k
      endif

      if(ntetype.gt.0) then
      	read(52,*) ttt,(ath(i),i=1,ntetype)
      	if(it.eq.iths+1.and.abs(ttt-time).gt.1.e-4) then
          write(11,*)'Starting time wrong for temp',it,ttt
          stop
      	endif
        
      	icount=0
      	do k=1,nope
          if(itetype(k).eq.1) then
            icount=icount+1
            if(icount.gt.ntetype) then
              write(11,*)'Wrong counting 3'
              stop
            endif
            tth(k)=ath(icount)
          endif
      	enddo !k
      endif

      if(nsatype.gt.0) then
      	read(53,*) ttt,(ath(i),i=1,nsatype)
      	if(it.eq.iths+1.and.abs(ttt-time).gt.1.e-4) then
          write(11,*)'Starting time wrong for salt',it,ttt
          stop
        endif

        icount=0
        do k=1,nope
          if(isatype(k).eq.1) then
            icount=icount+1
            if(icount.gt.nsatype) then
              write(11,*)'Wrong counting 4'
              stop
            endif
            sth(k)=ath(icount)
          endif
      	enddo !k
      endif

!...  Compute new vel. for flow b.c.
!...
      do i=1,nope
        if(ifltype(i)/=0) then
          isd1=iosd(i,1)
          n1=isidenode(isd1,1)
          n2=isidenode(isd1,2)
          H2=dps(isd1)+(peta(n1)+peta(n2))/2
          if(H2.lt.h0) then
            write(11,*)'Dry bnd side',n1,n2,H2
            stop
          endif
          uth(i)=qth(i)*ramp/cwidth(i)/H2*snx(isd1)
          vth(i)=qth(i)*ramp/cwidth(i)/H2*sny(isd1)
      	endif !ifltype(i)
      enddo !i=1,nope

      if(nscreen.eq.1) write(*,*)'done flow b.c.'
      write(16,*)'done flow b.c.'

!
!************************************************************************
!									*
!			Backtracking 					*
!									*
!************************************************************************
!

!...  Compute # of subdivisions for nsubfl=2
!...
      if(nsubfl==2) then
      	do i=1,np
          do k=1,nvrt
            sum=0
            do j=1,nne(i)
              ie=ine(i,j)
	      id=iself(i,j)
	      devm=0
	      do l=1,2
		if(l==1) then
	          nd=nm(ie,nx(i34(ie),id,1))
		else
	          nd=nm(ie,nx(i34(ie),id,i34(ie)-1))
		endif
	        rl=dsqrt((x(i)-x(nd))**2+(y(i)-y(nd))**2)
	        uvel1=(uu2(i,k)+uu2(i,k-1))/2
	        uvel2=(uu2(nd,k)+uu2(nd,k-1))/2
	        vvel1=(vv2(i,k)+vv2(i,k-1))/2
	        vvel2=(vv2(nd,k)+vv2(nd,k-1))/2
                dev=dsqrt((uvel1-uvel2)**2+(vvel1-vvel2)**2)/rl*dt
		if(dev>devm) devm=dev
	      enddo !l
              wfc=(ww2(i,k)+ww2(i,k-1))/2
              iw=2*wfc*dt/delzmin
              idelta=max0(int(devm/1.e-1),iw)
              sum=sum+dble(idelta)/nne(i)
            enddo !j=1,nne(i)
            ndelt(i,k)=max0(ndelt_min,min0(ndelt_max,int(sum)))
          enddo !k=1,nvrt
      	enddo !i=1,np
      endif !nsubfl=2

      st=0 !secnds(0.0) !timing the process

!...  From centers of faces 
!...
      do 451 i=1,ns
        n1=isidenode(i,1)
      	n2=isidenode(i,2)
      	if(kfe(is(i,1)).ne.0) then
          ie0=is(i,1)
      	else if(is(i,2).ne.0.and.kfe(is(i,2)).ne.0) then
          ie0=is(i,2)
      	else !do not btrack
          do j=kbs(i),kfs(i)
            vnbt(i,j)=vn2(i,j)
            vtbt(i,j)=vt2(i,j)
          enddo !j
          go to 451
        endif

        do j=kbs(i),kfs(i) !wet side
!         Initialize (x0,y0,z0),nnel and vel.
!	  Caution! nnel must be initialized inside this loop as it is updated inside.
          nnel=ie0
          jlev=j
          x0=xcj(i)
          y0=ycj(i)
          if(j.eq.kfe(nnel)) then
            z0=zmsl+eta1(nnel)-dc(nnel,j)/2
          else
            z0=ztot(j)-dc(nnel,j)/2
          endif
!	Scopes for ww2 etc. have been extended
          uuint=vn2(i,j)*snx(i)-vt2(i,j)*sny(i)
          vvint=vn2(i,j)*sny(i)+vt2(i,j)*snx(i)
          wdown=(ww2(n1,j-1)+ww2(n2,j-1))/2
          wup=(ww2(n1,j)+ww2(n2,j))/2
          wwint=(wdown+wup)/2
          vmag=dsqrt(uuint**2+vvint**2+wwint**2)

          if(vmag.le.1.e-4) then
!	    No activity
            vnbt(i,j)=vn2(i,j)
            vtbt(i,j)=vt2(i,j)
!	    Bottom friction
            if(j.eq.kbs(i)) tk(i)=Cd(i)*dsqrt(uuint**2+vvint**2)
          else !do backtrack
            ndels=(ndelt(n1,j)+ndelt(n2,j))/2
      	    dtb=dt/ndels !sub-step in backtracking
            call btrack(1,ndels,dtb,uuint,vvint,wwint,&
     &x0,y0,z0,xt,yt,zt,nnel,jlev,ttint,ssint,i)

            vnbt(i,j)=uuint*snx(i)+vvint*sny(i)
            vtbt(i,j)=-uuint*sny(i)+vvint*snx(i)
!         Bottom friction
            if(j.eq.kbs(i)) tk(i)=Cd(i)*dsqrt(uuint**2+vvint**2)
          endif

!         Overwrite if advection here is turned off
          if(nadv.eq.0.and.dabs(advt(i)).lt.0.1) then !advt=0
      	    vnbt(i,j)=vn2(i,j)
      	    vtbt(i,j)=vt2(i,j)
          endif
        enddo !j=kbs(i),kfs(i)
451   continue !i=1,ns

      en=0 !secnds(st)
      if(nscreen.eq.1) write(*,*)'backtracking took',en,'seconds...'
      write(16,*)'backtracking took',en,'seconds...'


!
!************************************************************************
!									*
!		Turbulence closure schemes				*
!	Compute turbulence diffusivities vdiff, tdiff,        		*
!	and in MY-G, also qdiff.					*
!									*
!************************************************************************
!

!... Richardson number at side and whole levels
      do i=1,ns
      	do k=kbs(i),kfs(i) !wet sides
          if(k.eq.kfs(i)) then 
            rich(i,k)=0
          else !k <= mj-1; implies M > m
            drhodz=(srho(i,k+1)-srho(i,k))/dzhalf(i,k)
            bvf=-g*(drhodz/rho0+g/1.5e3**2)
	    dundz=(vn2(i,k+1)-vn2(i,k))/dzhalf(i,k)
	    dutdz=(vt2(i,k+1)-vt2(i,k))/dzhalf(i,k)
            shear2=dundz**2+dutdz**2 !same form in 2 frames
            shear2=dmax1(shear2,1.0d-10)
            rich(i,k)=dmax1(bvf/shear2,0.0d0)
          endif
 	enddo !k=kbs(i),kfs(i)	
      enddo !i=1,ns

!... Scheme 1: Step function
      if(itur.eq.1) then
        richc=0.25  !critical Richardson #
        vdiffc=1.0d-3 !sub-critical mixing coefficient
        tdiffc=1.0*vdiffc
        do i=1,ns
          do k=kbs(i),kfs(i) !wet sides
   	    if(rich(i,k).lt.richc) then
              vdiff(i,k)=vdiffc
              tdiff(i,k)=tdiffc
            else
              vdiff(i,k)=1.0d-5
              tdiff(i,k)=1.0d-5
            endif
          enddo !k
        enddo !i=1,ns

        if(nscreen.eq.1) write(*,*) 'done turbulence closure (step)...'
        write(16,*) 'done turbulence closure (step)...'
      endif !itur=1

!... Scheme 2: Pacanowski and Philander (1981)
      if(itur.eq.2) then
        do i=1,ns
	  if(dps(i).le.h1_pp) then
	    vmax=vdmax_pp1
	    vmin=vdmin_pp1
	  else if(dps(i).lt.h2_pp) then
	    vmax=vdmax_pp1+(vdmax_pp2-vdmax_pp1)*(dps(i)-h1_pp)/(h2_pp-h1_pp)
	    vmin=vdmin_pp1+(vdmin_pp2-vdmin_pp1)*(dps(i)-h1_pp)/(h2_pp-h1_pp)
	  else !dps >= h2
	    vmax=vdmax_pp2
	    vmin=vdmin_pp2
	  endif

          do k=kbs(i),kfs(i) !wet sides
!	    vmax >= vmin
            vdiff(i,k)=vmax/(1+5*rich(i,k))**2+vmin
            tdiff(i,k)=vdiff(i,k)/(1+5*rich(i,k))+tdmin_pp
 	  enddo !k	
        enddo !i=1,ns

        if(nscreen.eq.1) write(*,*) 'done turbulence closure (PP)...'
        write(16,*) 'done turbulence closure (PP)...'
      endif !itur=2

!... Scheme 3: Mellor-Yamada-Galperin & Umlauf-Burchard scheme
      if(itur==3) then
!------------------------------------------------------------
      st=0 !secnds(0.0) !timing the process
!...
!... Solve for q2 and q2l (or psi in UB)
!...
!     Compute explicit terms
      do j=1,ns
        if(kfs(j)==0) cycle
        n1=isidenode(j,1)
        n2=isidenode(j,2)
        do k=kbs(j),kfs(j)
          q2bt(j,k)=(q2(j,k-1)+q2(j,k))/2
          xlbt(j,k)=(xl(j,k-1)+xl(j,k))/2
          if(xlbt(j,k).lt.xlmin2(j)) then
            write(11,*)'xlbt < xlmin2',it,j,k,xlbt(j,k),xlmin2(j)
            do l=kbs(j)-1,kfs(j)
              write(11,*)l,xlbt(j,l)
            enddo 
            stop
          endif
          udown=(uu2(n1,k-1)+uu2(n2,k-1))/2
          uup=(uu2(n1,k)+uu2(n2,k))/2
          vdown=(vv2(n1,k-1)+vv2(n2,k-1))/2
          vup=(vv2(n1,k)+vv2(n2,k))/2
          dudz=(uup-udown)/dz(j,k)
          dvdz=(vup-vdown)/dz(j,k)
          shearbt(j,k)=dudz**2+dvdz**2
          if(k==kfs(j).or.k==kbs(j)) then
            drhodz=0
          else ! m+1 <= k <= M-1
            rho_up=srho(j,k)+dz(j,k)/2/dzhalf(j,k)*(srho(j,k+1)-srho(j,k))
            rho_down=srho(j,k-1)+dz(j,k-1)/2/dzhalf(j,k-1)*(srho(j,k)-srho(j,k-1))
     	    drhodz=(rho_up-rho_down)/dz(j,k)
          endif
          rzbt(j,k)=g/rho0*drhodz
        enddo !k=kbs(j),kfs(j)
      enddo !j=1,ns


!	if(it.eq.1739) write(89,*)'WOW0',it

      do j=1,ns
      	if(kfs(j)==0) cycle

!	Compute upper bound for xl (not used)
      	if(is(j,2)/=0) then
          etam=dmax1(eta1(is(j,1)),eta1(is(j,2)))
      	else
          etam=eta1(is(j,1))
      	endif

      	do k=kbs(j),kfs(j)
          if(k==kfs(j)) then
            zctr=zmsl+etam-dz(j,k)/2
          else !M > m
            zctr=ztot(k)-dz(j,k)/2
          endif
          dists=zmsl+etam-zctr
          distb=zctr-(zmsl-dps(j))
!	  xlmax(k)=0.4*dmin1(dists,distb)
          xlmax(k)=0.4*dists*distb/(dps(j)+etam)
          if(xlmax(k)<=0) then
            write(11,*)'Dist<0 in MY-G',j,k,dps(j)+etam,dists,distb
            stop
          endif
        enddo !k

!	if(it.eq.1739) write(89,*)'WOW1',it,j

!	b.c. (computed using values from previous time except wind)
        q2fs=16.6**(2.0/3)*dsqrt(tau_n(j)**2+tau_t(j)**2)/2
      	q2fs=dmax1(q2fs,q2min)
      	q2bot=16.6**(2.0/3)*tk(j)*dsqrt(vn2(j,kbs(j))**2+vt2(j,kbs(j))**2)/2
      	q2bot=dmax1(q2bot,q2min)
        xlfs=dmax1(xlmin1(j),dz(j,kfs(j))/2*0.4)
      	xlbot=dmax1(xlmin2(j),dz(j,kbs(j))/2*0.4)

      	nqdim=kfs(j)-kbs(j)+1 !>=1

!	if(it.eq.1739) write(89,*)'WOW2',it,j

!...    Need to solve eq. system since nqdim >=1 (or M >= m)
!	Matrix Q
      	do k=1,nqdim !row #
          klev=kfs(j)-k+1 !level #; m <= klev <= M
	  eleft=0 !for bdia
	  eright=0
          if(k>1) then
            eleft=-qdiff(j,klev)*dt/dzhalf(j,klev)
	    alow(k)=eleft
	  endif
          if(k<nqdim) then
            eright=-qdiff(j,klev-1)*dt/dzhalf(j,klev-1)
	    cupp(k)=eright
	  endif

	  if(k==nqdim) then
	    prod=vdiff(j,klev)*shearbt(j,klev)
	    buoy=tdiff(j,klev)*rzbt(j,klev)
	  else
	    prod=(vdiff(j,klev)+vdiff(j,klev-1))/2*shearbt(j,klev)
	    buoy=(tdiff(j,klev)+tdiff(j,klev-1))/2*rzbt(j,klev)
	  endif
	  diss=cmiu0**3*dsqrt(q2bt(j,klev))/xlbt(j,klev) !diss/k
	  if(prod+buoy>0) then
	    pplus=prod+buoy
	    pminus=0
	  else
	    pplus=0
	    pminus=-prod-buoy
	  endif
          bdia(k)=dz(j,klev)-eleft-eright+(diss+pminus/q2bt(j,klev))*dt*dz(j,klev)

!	  R.H.S. for q2
    	  rrhs(k,1)=dz(j,klev)*(q2bt(j,klev)+dt*pplus)
        enddo !k=1,nqdim

!	Soln for q2 at new level
      	call tridag(mnv,nqdim,1,alow,bdia,cupp,rrhs,soln,gam)
        do k=1,nqdim
          klev=kfs(j)-k+1 !level #; m <= klev <= M
	  if(klev==kfs(j)) then !fs value previals for M=m
      	    q2tmp(klev)=q2fs
	  else if(klev==kbs(j)) then
      	    q2tmp(klev)=q2bot
	  else
            q2tmp(klev)=dmax1(soln(k,1),q2min)
	  endif
        enddo !k

!	if(it.eq.1739) write(89,*)'WOW4',it,j

!	Matrix QL
      	do k=1,nqdim !row #
          klev=kfs(j)-k+1 !level #; m <= klev <= M
	  eleft=0
	  eright=0
          if(k>1) then !klev <= M-1
            eleft=-qdiff2(j,klev)*dt/dzhalf(j,klev)
	    alow(k)=eleft
	  endif
          if(k<nqdim) then !klev >= m+1
            eright=-qdiff2(j,klev-1)*dt/dzhalf(j,klev-1)
	    cupp(k)=eright
	  endif

!	  Diagonal and RHS
	  bdia(k)=dz(j,klev)-eleft-eright
	  if(mid.eq.'MY'.or.mid.eq.'KL') then
            zctr=ztot(klev)-dz(j,klev)/2
            dists=dmax1(zmsl+etam-zctr,xlmin2(j))
            distb=dmax1(zctr-(zmsl-dps(j)),xlmin2(j))
            psi=1+1.33*(xlbt(j,klev)/0.4/distb)**2+0.25*(xlbt(j,klev)/0.4/dists)**2
	    cpsi2p=psi*cpsi2 !F_wall*cpsi2
	  else !other GLS
	    cpsi2p=cpsi2
	  endif
	  diss=cpsi2p*cmiu0**3*dsqrt(q2bt(j,klev))/xlbt(j,klev) !diss/k
	  bdia(k)=bdia(k)+dt*dz(j,klev)*diss

	  if(k==nqdim) then
	    prod=cpsi1*vdiff(j,klev)*shearbt(j,klev)
	  else
	    prod=cpsi1*(vdiff(j,klev)+vdiff(j,klev-1))/2*shearbt(j,klev)
	  endif
	  if(mid.eq.'MY') then
	    cpsi3=0.9
	  else !GLS models
	    if(rzbt(j,klev)>0) then !unstable
	      cpsi3=1
	    else !stable
	      select case(mid)
	        case('KL')
		  cpsi3=2.53
	        case('KE')
		  cpsi3=-0.52
	        case('KW')
		  cpsi3=-0.58
	        case('UB')
		  cpsi3=0.1
	        case default
	          write(11,*)'Unknown closure model:',mid
	  	  stop
	      end select
	    endif
	  endif

	  if(k==nqdim) then
	    buoy=cpsi3*tdiff(j,klev)*rzbt(j,klev)
	  else
	    buoy=cpsi3*(tdiff(j,klev)+tdiff(j,klev-1))/2*rzbt(j,klev)
	  endif
	  if(prod+buoy>0) then
	    pplus=prod+buoy
	    pminus=0
	  else
	    pplus=0
	    pminus=-prod-buoy
	  endif
          bdia(k)=bdia(k)+pminus/q2bt(j,klev)*dt*dz(j,klev)
          if(klev==kbs(j).or.klev==kfs(j)) bdia(k)=bdia(k)+dt*0.4*rnub*qdiff2(j,klev)/xlbt(j,klev)

	  tmp=cmiu0**rpub*q2bt(j,klev)**rmub*xlbt(j,klev)**rnub !psi^n_{j,k}
	  rrhs(k,1)=dz(j,klev)*tmp*(1+dt*pplus/q2bt(j,klev))
      	enddo !k=1,nqdim

!	Soln for q2l and xl at new level
        call tridag(mnv,nqdim,1,alow,bdia,cupp,rrhs,soln,gam)
        do k=1,nqdim
          klev=kfs(j)-k+1
          q2l=dmax1(soln(k,1),psimin)
	  if(klev==kfs(j)) then
            xltmp(klev)=xlfs
	  else if(klev==kbs(j)) then
            xltmp(klev)=xlbot
	  else
            xltmp(klev)=(q2l*cmiu0**(-rpub)*q2tmp(klev)**(-rmub))**(1/rnub)
	  endif
!	  Galperin's clipping 
	  if(rzbt(j,klev)<0) then
	    upper=dsqrt(-0.56*q2tmp(klev)/rzbt(j,klev))
	    xltmp(klev)=dmin1(xltmp(klev),upper)
	  endif
!	  Max. length based on dissipation; xlmin2 prevails
	  xl_max=(cmiu0*dsqrt(q2tmp(klev)))**3/eps_min
          xltmp(klev)=dmax1(xlmin2(j),dmin1(xl_max,xltmp(klev)))
        enddo !k

!	if(it.eq.1739) write(89,*)'WOW6',it,j

!	Convert to whole levels
        do k=kbs(j)-1,kfs(j)
          if(k==kbs(j)-1) then
            q2(j,k)=q2bot
            xl(j,k)=xlbot
          else if(k==kfs(j)) then
            q2(j,k)=q2fs
            xl(j,k)=xlfs
          else !M > m and m<= k <= M-1
            q2(j,k)=q2tmp(k)+dz(j,k)/2/dzhalf(j,k)*(q2tmp(k+1)-q2tmp(k))
            xl(j,k)=xltmp(k)+dz(j,k)/2/dzhalf(j,k)*(xltmp(k+1)-xltmp(k))
            if(q2(j,k)<0.or.xl(j,k)<xlmin2(j)) then
              write(11,*)'Negative q2/xl',q2(j,k),xl(j,k)
              stop
            endif
          endif
      	enddo !k=kbs(j)-1,kfs(j)

      enddo !j=1,ns

!      if(it.eq.1739) write(89,*)'WOW7',it

!...
!... Compute vertical diffusivities at new time
!...
      do i=1,ns
      	do j=kbs(i),kfs(i) !wet sides
	  call asm(g,i,j,vd,td,qd,qd2)
          vdiff(i,j)=dmin1(diffmax(i),dmax1(bgdiff,vd))
          tdiff(i,j)=dmin1(diffmax(i),dmax1(bgdiff,td))
          qdiff(i,j)=dmin1(diffmax(i),dmax1(bgdiff,qd))
          qdiff2(i,j)=dmin1(diffmax(i),dmax1(bgdiff,qd2))
        enddo !j=kbs(i),kfs(i)
      enddo !i=1,ns

!...  Extend q2 and xl etc. to account for level changes next time
      do i=1,ns
      	if(kfs(i).eq.0) cycle
      	do j=0,kbs(i)-2
          q2(i,j)=q2(i,kbs(i)-1)
          xl(i,j)=xl(i,kbs(i)-1)
      	enddo 
      	do j=1,kbs(i)-1
          tdiff(i,j)=tdiff(i,kbs(i))
          vdiff(i,j)=vdiff(i,kbs(i))
          qdiff(i,j)=qdiff(i,kbs(i))
          qdiff2(i,j)=qdiff2(i,kbs(i))
        enddo !j=1,kbs(i)-1
        do j=kfs(i)+1,nvrt
          q2(i,j)=q2(i,kfs(i))
          xl(i,j)=xl(i,kfs(i))
          tdiff(i,j)=tdiff(i,kfs(i))
          vdiff(i,j)=vdiff(i,kfs(i))
          qdiff(i,j)=qdiff(i,kfs(i))
          qdiff2(i,j)=qdiff2(i,kfs(i))
        enddo !j=kfs(i)+1,nvrt
      enddo !i

      en=0 !secnds(st)
      if(nscreen.eq.1) write(*,*)'MYG-UB took',en,'seconds'
      write(16,*)'MYG-UB took',en,'seconds'
!------------------------------------------------------------
      endif !itur=3
!...
!... End of MY-G scheme

!... compute diffusivity at point and whole level: tdiffp for transport eq.
      do i=1,np
        do k=kbp(i),kfp(i) !wet nodes only
          tdiffp(i,k)=0
          icount=0
          do j=1,nne(i)
            ie=ine(i,j)
            do l=1,i34(ie)
              isd=js(ie,l)
              if(kfs(isd)/=0) then
                icount=icount+1
                klev=min(max(k,kbs(isd)),kfs(isd))
                tdiffp(i,k)=tdiffp(i,k)+tdiff(isd,klev)
              endif
            enddo !l
          enddo !j

          if(icount==0) then
            tdiffp(i,k)=1.e-2
          else
            tdiffp(i,k)=tdiffp(i,k)/icount
          endif
        enddo !k=kbp(i),kfp(i)
!...	Extend
	do k=kfp(i)+1,nvrt
	  tdiffp(i,k)=tdiffp(i,kfp(i))
	enddo !k
      enddo !i=1,np

!... End turbulence closure schemes


!
!************************************************************************
!									*
!		Momentum equations					*
!									*
!************************************************************************
!

      st=0 !secnds(0.0) !timing the process
!... compute matrices and vectors for momentum equation, along each side
!... zhat1(*,*) is Ainv*DeltaZ
!... ghat(*,*) is Ainv*gvec
!... gvec(*) is the vector G
!... 
      do 419 j=1,ns
      	if(kfs(j).eq.0) go to 419

      	node1=isidenode(j,1)
      	node2=isidenode(j,2)
      	nv=kfs(j)-kbs(j)+1

!	Coefficient matrix
!	Middle rows
        do k=kbs(j)+1,kfs(j)-1
          itmp=kfs(j)-k+1 !row #
          alow(itmp)=-vdiff(j,k)*dt/dzhalf(j,k) 
          cupp(itmp)=-vdiff(j,k-1)*dt/dzhalf(j,k-1)
          bdia(itmp)=dz(j,k)-alow(itmp)-cupp(itmp)
         enddo !k

!	1st and last rows
        if(kfs(j).eq.kbs(j)) then
          bdia(1)=dz(j,kfs(j))+tk(j)*dt
        else !M > m
          alow(nv)=-vdiff(j,kbs(j))*dt/dzhalf(j,kbs(j))
          bdia(nv)=dz(j,kbs(j))-alow(nv)+tk(j)*dt
          cupp(1)=-vdiff(j,kfs(j)-1)*dt/dzhalf(j,kfs(j)-1)
          bdia(1)=dz(j,kfs(j))-cupp(1)
        endif

!	Pre-compute tidal potential gradients for each side to save time
	  if(ntip.gt.0.and.dps(j).ge.tip_dp.and.is(j,2).ne.0) then
	    do l=1,2
	      ie=is(j,l)
	      swild(l)=0
	      do jf=1,ntip
	 	ncyc=int(tfreq(jf)*time/2/pi)
              	arg=tfreq(jf)*time-ncyc*2*pi+jspc(jf)*xlon_e(ie)+tear(jf)
              	swild(l)=swild(l)+ramptide*tamp(jf)*tnf(jf)*fun_lat_e(ie,jspc(jf))*dcos(arg)	
	      enddo !jf
	    enddo !l=1,2
	    detp_dx=(swild(2)-swild(1))/delj(j)
	  else
	    detp_dx=0
	  endif !ntip.gt.0.and.is(j,2).ne.0

	  if(ntip.gt.0.and.dps(j).ge.tip_dp) then
	    do l=1,2
	      nd=isidenode(j,l)
	      swild(l)=0
	      do jf=1,ntip
	 	ncyc=int(tfreq(jf)*time/2/pi)
              	arg=tfreq(jf)*time-ncyc*2*pi+jspc(jf)*xlon(nd)+tear(jf)
              	swild(l)=swild(l)+ramptide*tamp(jf)*tnf(jf)*fun_lat(nd,jspc(jf))*dcos(arg)	
	      enddo !jf
	    enddo !l=1,2
	    detp_dy=(swild(2)-swild(1))/distj(j)
	  else
	    detp_dy=0
	  endif !ntip.gt.0

!	Vectors F & G
        do k=1,kfs(j)-kbs(j)+1
          klev=kfs(j)+1-k

!	  b.c. to be imposed at the end

!	  Coriolis, advection and wind stress
          gvec(k)=dz(j,klev)*(vnbt(j,klev)+cori(j)*dt*vt2(j,klev))
          fvec(k)=dz(j,klev)*(vtbt(j,klev)-cori(j)*dt*vn2(j,klev))
          if(k.eq.1) then
            gvec(k)=gvec(k)+dt*tau_n(j)
            fvec(k)=fvec(k)+dt*tau_t(j)
          endif

!	  Elevation gradient
          if(is(j,2).ne.0) gvec(k)=gvec(k)-dz(j,klev)*g*dt*(1-thetai)* &
     &(eta1(is(j,2))-eta1(is(j,1)))/delj(j)
          if(ibad(node1).eq.0.and.ibad(node2).eq.0) fvec(k)=fvec(k)- &
     &dz(j,klev)*g*dt*(1-thetai)*(peta(node2)-peta(node1))/distj(j)

!	  Earth tidal potential
	  gvec(k)=gvec(k)+0.69*dz(j,klev)*g*dt*detp_dx
	  fvec(k)=fvec(k)+0.69*dz(j,klev)*g*dt*detp_dy

!	  Atmos. pressure gradient
          fvec(k)=fvec(k)-dz(j,klev)*dt/rho0*(pr1(node2)-pr1(node1))/distj(j)
          if(is(j,2).ne.0) then
            do l=1,2
	      swild(l)=0
              iel=is(j,l)
	      do ll=1,i34(iel)
	        swild(l)=swild(l)+pr1(nm(iel,ll))/i34(iel)
	      enddo !ll
            enddo !l=1,2
            gvec(k)=gvec(k)-dz(j,klev)*dt/rho0*(swild(2)-swild(1))/delj(j)
          endif

!         Horizontal diffusion
!...      Smagorinsky diffusivity at half levels 
!...
      	  if(smagcoef.ne.0.and.is(j,2).ne.0) then
	    do l=1,2
              nd=isidenode(j,l)
              utmp=(uu2(nd,klev)+uu2(nd,klev-1))/2
              vtmp=(vv2(nd,klev)+vv2(nd,klev-1))/2
	      swild(2*l-1)=utmp*snx(j)+vtmp*sny(j) !uvel_n[12]l
	      swild(2*l)=-utmp*sny(j)+vtmp*snx(j) !vvel_n[12]l
	    enddo !l

	    do l=1,2
              ie=is(j,l)
	      utmp=0
	      vtmp=0
	      do ll=1,i34(ie)
	        nd=nm(ie,ll)
	        utmp=utmp+(uu2(nd,klev)+uu2(nd,klev-1))/i34(ie)/2
	        vtmp=vtmp+(vv2(nd,klev)+vv2(nd,klev-1))/i34(ie)/2
	      enddo !ll
	      swild(2*l+3)=utmp*snx(j)+vtmp*sny(j) !uvel_e[12]l
	      swild(2*l+4)=-utmp*sny(j)+vtmp*snx(j) !vvel_e[12]l
	    enddo !l

            if(ihorcon.eq.0) then
              hdiff=smagcoef
            else !Smagorinsky
              dudx=(swild(7)-swild(5))/delj(j)
              dvdx=(swild(8)-swild(6))/delj(j)
              dudy=(swild(3)-swild(1))/distj(j)
              dvdy=(swild(4)-swild(2))/distj(j)
              hdiff=smagcoef*areas(is(j,1))*dsqrt(dudx**2+dvdy**2+(dvdx+dudy)**2/2)
            endif

      	    d2udy=4*(swild(3)+swild(1)-2*vn2(j,klev))/distj(j)**2
      	    d2vdy=4*(swild(4)+swild(2)-2*vt2(j,klev))/distj(j)**2

	    do l=1,2
	      ie=is(j,l)
              rl=dsqrt((xctr(ie)-xcj(j))**2+(yctr(ie)-ycj(j))**2)
              swild(2*l+7)=(swild(2*l+3)-vn2(j,klev))/rl
              swild(2*l+8)=(swild(2*l+4)-vt2(j,klev))/rl
	    enddo !l
            d2udx=2*(swild(9)+swild(11))/delj(j)
            d2vdx=2*(swild(10)+swild(12))/delj(j)

            gvec(k)=gvec(k)+dz(j,klev)*dt*hdiff*(d2udx+d2udy)
            fvec(k)=fvec(k)+dz(j,klev)*dt*hdiff*(d2vdx+d2vdy)
          endif !is(j,2).ne.0

!	  Baroclinic gradient term 
          if(ibc.eq.0) then
            do l=klev,kfs(j)
	      if(is(j,2).ne.0.and.l.ge.kbe(is(j,1)).and.l.le.kfe(is(j,1)).and.l.ge.kbe(is(j,2)).&
     &and.l.le.kfe(is(j,2))) then
	        tmp=-dz(j,klev)*g*dt/rho0*dz(j,l)*(erho(is(j,2),l)-erho(is(j,1),l))/delj(j)
                if(l.eq.klev) then
                  gvec(k)=gvec(k)+tmp/2*rampbc
                else
                  gvec(k)=gvec(k)+tmp*rampbc
                endif
	      endif

	      if(l.ge.kbp(node1).and.l.le.kfp(node1).and.&
     &l.ge.kbp(node2).and.l.le.kfp(node2)) then !2 nodes wet
	        tmp=-dz(j,klev)*g*dt/rho0*dz(j,l)*(prho(node2,l)-prho(node1,l))/distj(j)
                if(l.eq.klev) then
                  fvec(k)=fvec(k)+tmp/2*rampbc
                else
                  fvec(k)=fvec(k)+tmp*rampbc
                endif
	      endif 
	    enddo !l=klev,kfs(j)
          endif !ibc.eq.0

!	  Impose b.c.
          if(isflowside(j).gt.0) then
!	    Impose flux b.c.
            ibnd=isflowside(j)
            gvec(k)=uth(ibnd)*snx(j)+vth(ibnd)*sny(j)
            fvec(k)=-uth(ibnd)*sny(j)+vth(ibnd)*snx(j)
     	  else if(isbs(j).eq.-1) then
!	    Land b.c.
            gvec(k)=0
          endif !b.c.
        enddo !k=1,nv

!	RHS
        do k=1,nv
          rrhs(k,1)=gvec(k)
          rrhs(k,2)=fvec(k)
          rrhs(k,3)=dz(j,kfs(j)+1-k)
        enddo !k

!	Causion: the eq. is wrong for b.c.!
      	call tridag(mnv,nv,3,alow,bdia,cupp,rrhs,soln,gam)

!	Fhat and Ghat
        do k=1,nv
          if(isflowside(j).gt.0.or.isbs(j).eq.-1) then
            ghat(j,k)=gvec(k)
          else
            ghat(j,k)=soln(k,1)
          endif

          if(isflowside(j).gt.0) then
            fhat(j,k)=fvec(k)
          else
            fhat(j,k)=soln(k,2)
          endif
        enddo !k=1,nv

!	zhat[1,2]
        do k=1,nv
          if(isflowside(j).gt.0.or.is(j,2).eq.0) then
            zhat1(j,k)=0
          else
            zhat1(j,k)=soln(k,3)
          endif

          if(isflowside(j).gt.0) then
            zhat2(j,k)=0
          else
            zhat2(j,k)=soln(k,3)
          endif
        enddo !k
419   continue !j=1,ns

      en=0 !secnds(st)
      if(nscreen.eq.1) write(*,*)'setting up momentum eq took',en,'sec'
      write(16,*) 'setting up momentum eq took',en,'sec'

!
!************************************************************************
!									*
!		Wave-continuity equations				*
!									*
!************************************************************************
!

!...  compute elevation boundary conditions
!...
      do i=1,nope
        do j=1,noe(i)
          iel=ioe(i,j)
          if(iettype(i).eq.1.or.iettype(i).eq.2) then
            elbc(iel)=ramp*eth(i)
          else if(iettype(i).eq.3) then
            elbc(iel)=0 !initialize
            do jfr=1,nbfr
	      ncyc=int(amig(jfr)*time/2/pi)
              arg=amig(jfr)*time-ncyc*2*pi+face(jfr)-efa(i,jfr,j)
              elbc(iel)=elbc(iel)+ramptide*ff(jfr)*emo(i,jfr,j)*dcos(arg)
            enddo !jfr=1,nbfr
          else if(iettype(i).eq.-1) then
            sum=0
            icount=0
            do k=1,i34(iel)
              nd=nm(iel,k)
              do l=1,nne(nd)
                ie=ine(nd,l)
                if(isbe(ie).eq.0) then
                  icount=icount+1
                  sum=sum+eta1(ie)
               endif
              enddo !l
            enddo !k

            if(icount.eq.0) then
              write(11,*)'Isolated obe cannot have obc',iel
              stop
            else
              elbc(iel)=sum/icount
            endif
          endif
        enddo !j=1,noe(i)
      enddo !i=1,nope

!...  setup coefficient matrix, sparsem, for the wave equation
!...  setup RHS, qel, for the wave equation using etmp1, and etmp3
!...  No b.c. are involved yet
      do i=1,ne
        sparsem(i,0)=areas(i)
        do j=1,i34(i)
          nel=ic3(i,j)
          jsj=js(i,j)
          sparsem(i,j)=0 !in case if below is not true
   	  if(kfs(jsj)/=0.and.nel/=0) then
!	    DZ*etmp1
            const=0
            do k=1,kfs(jsj)-kbs(jsj)+1
              const=const+dz(jsj,kfs(jsj)+1-k)*zhat1(jsj,k)
            enddo
            const=g*dt**2*thetai**2*distj(jsj)/delj(jsj)*const  
            if(const.lt.0) then
              write(11,*)'Not positive definite',i,const
              stop
            endif
            sparsem(i,0)=sparsem(i,0)+const
            sparsem(i,j)=-const
	  else if(isflowside2(jsj)>0) then
!	  Flather's b.c.
	    const=0
	    do k=kbs(jsj),kfs(jsj)
	       const=const+dt*distj(jsj)*dsqrt(g/dps(jsj))*dz(jsj,k)
	    enddo !k
	    sparsem(i,0)=sparsem(i,0)+const
          endif 
        enddo !j=1,i34(i)

!	qel
        qel(i)=eta1(i)*areas(i)
        do j=1,i34(i)
          jsj=js(i,j)
          if(kfs(jsj)/=0) then
            const1=0
            const2=0
            do k=1,kfs(jsj)-kbs(jsj)+1
              klev=kfs(jsj)+1-k
              const1=const1+dz(jsj,klev)*vn2(jsj,klev)
              const2=const2+dz(jsj,klev)*ghat(jsj,k) 
            enddo
            qel(i)=qel(i)-(1-thetai)*dt*ssign(i,j)*distj(jsj)*const1- &
     &thetai*dt*ssign(i,j)*distj(jsj)*const2
          endif

!	  Flather's b.c.
	  if(isflowside2(jsj)>0) then
	    ibnd=isflowside2(jsj)
	    vnorm=uth(ibnd)*snx(jsj)+vth(ibnd)*sny(jsj)
	    do k=kbs(jsj),kfs(jsj)
	      qel(i)=qel(i)-dt*distj(jsj)*dz(jsj,k)*vnorm
	    enddo !k
	  endif
        enddo !j=1,i34
      enddo !i=1,ne

!...  Create a mapping between element index and actual eq. index
      icount=0
      do i=1,ne
        if(isbe(i)/=0.and.iettype(isbe(i))/=0) then
          icount=icount+1
        else
          imape(i)=i-icount !only for non-bnd elements
        endif
      enddo !i

!...  Assemble sparse matrix format
!...
      neq=0 !final index of eqs.
      nnz=0 !!# of non-zero entries
      do i=1,ne
        if(isbe(i)/=0.and.iettype(isbe(i))/=0) cycle

        neq=neq+1
        nnz=nnz+1
        icoef(neq)=nnz
        jcoef(nnz)=imape(i)
        if(imape(i)/=neq) then
          write(11,*)'Impossible 300',i
          stop
        endif
        e2coef(nnz)=sparsem(i,0)
        qel2(neq)=qel(i)
        eta3(neq)=eta2(i) !initial guess
        do j=1,i34(i)
          nel=ic3(i,j)
          if(nel/=0) then
            if(isbe(nel)/=0.and.iettype(isbe(nel))/=0) then !b.c.
              qel2(neq)=qel2(neq)-sparsem(i,j)*elbc(nel)
            else if(nel>i) then !upper triangle only
              nnz=nnz+1
              jcoef(nnz)=imape(nel)
              e2coef(nnz)=sparsem(i,j)
            endif
          endif
        enddo !j=1,i34(i)
      enddo !i=1,ne

      if(neq.le.0) then
        write(11,*)'No eqs. to be solved!'
        stop
      endif
      icoef(neq+1)=nnz+1


!...  solve the wave equation for elevations at each element center
!...  jcg jacobi conjugate gradient solver from itpack2d, srcv2d.f
!...
      st=0 !secnds(0.0) !timing the process

!...  input information about solver
!...
      call dfault(iparm,rparm)
      iparm(1)=itmax1
      iparm(2)=1 !level of output msg
      iparm(4)=33 !output msg to fort.??
      iparm(5)=0 !symmetric system
!      iparm(10)=iremove
      iparm(11)=1 !no timing
      iparm(12)=1 !error analysis
      rparm(1)=zeta  !0.11102230E-13 !stopping criterion
      rparm(8)=tol  !1.0e-8

!jcg    call nspcg (jac1,basic,ndim,mdim,nor,4,
!jcg &        e2coef,jcoef,p,ip,eta2,
!jcg &        ubar,qel,wksp,iwksp,nwksp,inw,iparm,rparm,ier)

      if(isolver.eq.1) then
        call jcg(neq,icoef,jcoef,e2coef,qel2,eta3,iwksp,nwksp,wksp,iparm,rparm,ier)
      else if(isolver.eq.2) then
        call jsi(neq,icoef,jcoef,e2coef,qel2,eta3,iwksp,nwksp,wksp,iparm,rparm,ier)
      else if(isolver.eq.3) then
        call ssorcg(neq,icoef,jcoef,e2coef,qel2,eta3,iwksp,nwksp,wksp,iparm,rparm,ier)
      else !isolver=4
        call ssorsi(neq,icoef,jcoef,e2coef,qel2,eta3,iwksp,nwksp,wksp,iparm,rparm,ier)
      endif

      en=0 !secnds(st)
      if(nscreen.eq.1) write(*,*)'calling solver took',en,'seconds...',iparm(1),iparm(8)
      write(16,*)'calling solver took',en,'seconds...',iparm(1),iparm(8)

!...  Re-assemble new elevations
!...
      do i=1,ne
        if(isbe(i)/=0.and.iettype(isbe(i))/=0) then
          eta2(i)=elbc(i)
      	else
          eta2(i)=eta3(imape(i))
      	endif
      enddo !i

!... Sponge layer
!...
      if(isponge.ne.0) then
        do i=1,ne
	  rel=0
	  do j=1,i34(i)
            rel=rel+relax(nm(i,j))/i34(i)
	  enddo !j
          eta2(i)=eta2(i)*rel+etaic(i)*(1-rel)
        enddo !i
      endif

!... Output total elevation
      etatot=0
      do i=1,ne
        etatot=etatot+dabs(eta2(i))
      enddo
      if(nscreen.eq.1) write(*,*) 'etatot=',etatot
      write(16,*) 'etatot=',etatot

!...  Compute nodal elevations for tang. vel.
      do i=1,np
        ibad(i)=0 !reliable elevation

        sum1=0 !sum of elemental elevations
        sum=0 !sum of areas
        do j=1,nne(i)
          ie=ine(i,j)
          if(kfe(ie).eq.0) then
            ibad(i)=1
          else
            sum1=sum1+areas(ie)*eta2(ie)
            sum=sum+areas(ie)
          endif
        enddo !j=1,nne(i)

        if(ibad(i).eq.0) peta(i)=sum1/sum
      enddo !i=1,np

!     Estimate bad elevations for rewetting
      do i=1,np
      	if(ibad(i).eq.1) then
          sum=0
          icount=0
          do j=1,nnp(i)
            nd=inp(i,j)
            if(ibad(nd).eq.0) then
              icount=icount+1
              sum=sum+peta(nd)
            endif
          enddo !j
          if(icount.ne.0) peta(i)=sum/icount
        endif !ibad(i).eq.1
      enddo !i=1,np

!...  solve the momentum equation for normal & tang. velocities 
!...
      do i=1,ns
        n1=isidenode(i,1)
        n2=isidenode(i,2)
        if(kfs(i).eq.0) then
          do k=0,nvrt
            vn2(i,k)=0
            vt2(i,k)=0
          enddo !k
        else !wet side
          do j=1,kfs(i)-kbs(i)+1
            jlev=kfs(i)+1-j
            vn2(i,jlev)=ghat(i,j)
       	    if(is(i,2).ne.0) vn2(i,jlev)=vn2(i,jlev)-thetai*g*dt/delj(i)* &
     &(eta2(is(i,2))-eta2(is(i,1)))*zhat1(i,j)
            vt2(i,jlev)=fhat(i,j)
            if(ibad(n1).eq.0.and.ibad(n2).eq.0) vt2(i,jlev)=vt2(i,jlev)-thetai*g*dt/distj(i)* &
     &(peta(n2)-peta(n1))*zhat2(i,j)
!	    Impose bounds for tang. vel. for initially dry nodes (dp < 0)
	    if(idrynode(n1).eq.1.or.idrynode(n2).eq.1) vt2(i,jlev)=dmax1(-5.d0,dmin1(vt2(i,jlev),5.d0))
          enddo !j=1,kfs(i)-kbs(i)+1

!	  Impose Flather's b.c.
	  if(isflowside2(i)>0) then
	    ibnd=isflowside2(i)
            vnorm=uth(ibnd)*snx(i)+vth(ibnd)*sny(i)
	    do k=kbs(i),kfs(i)
	      vt2(i,k)=0
	      vn2(i,k)=vnorm+dsqrt(g/dps(i))*eta2(is(i,1))
	    enddo !k
	  endif

!     	  Extend scope 
 	  do j=0,kbs(i)-1
    	    vn2(i,j)=0  !vn2(i,kbs(i))
    	    vt2(i,j)=0 
          enddo
          do j=kfs(i)+1,nvrt
            vn2(i,j)=vn2(i,kfs(i))
            vt2(i,j)=vt2(i,kfs(i))
          enddo
        endif !kfs(i).ne.0

      enddo !i=1,ns

      if(nscreen.eq.1) write(*,*)'done solving momentum eq...'
      write(16,*)'done solving momentum eq...'

!...  solve for vertical velocities
!...
      do 437 i=1,ne
 	if(kfe(i).eq.0) go to 437

 	we(i,kbe(i)-1)=0
        do l=kbe(i),kfe(i)
          we(i,l)=we(i,l-1)
          do j=1,i34(i)
            jsj=js(i,j)
            if(l.ge.kbs(jsj).and.l.le.kfs(jsj)) we(i,l)=we(i,l)-ssign(i,j)*distj(jsj)*&
     &dz(jsj,l)*vn2(jsj,l)/areas(i) !wet
          enddo !j=1,3
        enddo !l

!	Extend
        do l=0,kbe(i)-2
          we(i,l)=0
        enddo
        do l=kfe(i)+1,nvrt
          we(i,l)=we(i,kfe(i))
        enddo
437   continue  !i=1,ne

      if(nscreen.eq.1) write(*,*)'done solving w...'
      write(16,*)'done solving w...'

!     Test backtracking alone
!      do i=1,np
!         peta(i)=0
!       enddo
!       do i=1,ne
! 	eta2(i)=0
!         do k=0,nvrt
! 	  we(i,k)=0
!         enddo
!       enddo
!       do i=1,ns
!         do k=0,nvrt
! 	  vn2(i,k)=0.1*snx(i) 
! 	  vt2(i,k)=-0.1*sny(i)
!         enddo
!       enddo

!...  Recompute levels
      call levels(iths,it,itur)
      if(nscreen.eq.1) write(*,*) 'done recomputing levels...'
      write(16,*) 'done recomputing levels...'

!...  Compute nodal vel. for output and next backtracking
      call nodalvel

!...  Update elevation
!...
      do i=1,ne
        eta1(i)=eta2(i)
      enddo

!
!************************************************************************
!									*
!           Mode splitting: backtracking for transport			*
!									*
!************************************************************************
!
      if((ibc.eq.0.or.ibtp.eq.1).and.mod(it,nsplit).eq.0) then
!----------------------------------------------------------------------
      st=0 !secnds(0.0) !timing the process

!...  From centers of faces 
!...
      do 478 i=1,ns
        n1=isidenode(i,1)
        n2=isidenode(i,2)
        if(kfe(is(i,1)).ne.0) then
          ie0=is(i,1)
        else if(is(i,2).ne.0.and.kfe(is(i,2)).ne.0) then
          ie0=is(i,2)
        else !do not btrack
          do j=kbs(i),kfs(i)
            tsdbt(i,j)=tsd(i,j)
            ssdbt(i,j)=ssd(i,j)
          enddo !j
          go to 478
        endif

        do j=kbs(i),kfs(i) !wet side
!         Initialize (x0,y0,z0),nnel and vel.
!	  Caution! nnel must be initialized inside this loop as it is updated inside.
          nnel=ie0
          jlev=j
          x0=xcj(i)
          y0=ycj(i)
          if(j.eq.kfe(nnel)) then
            z0=zmsl+eta1(nnel)-dc(nnel,j)/2
          else
            z0=ztot(j)-dc(nnel,j)/2
          endif
!	Scopes for ww2 etc. have been extended
          uuint=vn2(i,j)*snx(i)-vt2(i,j)*sny(i)
          vvint=vn2(i,j)*sny(i)+vt2(i,j)*snx(i)
          wdown=(ww2(n1,j-1)+ww2(n2,j-1))/2
          wup=(ww2(n1,j)+ww2(n2,j))/2
          wwint=(wdown+wup)/2
          vmag=dsqrt(uuint**2+vvint**2+wwint**2)

          if(vmag.le.1.e-4) then
!	  No activity
            tsdbt(i,j)=tsd(i,j)
            ssdbt(i,j)=ssd(i,j)
          else !do backtrack
            ndels=(ndelt(n1,j)+ndelt(n2,j))/2
      	    dtb=dt/ndels !sub-step in backtracking
            call btrack(2,nsplit*ndels,dtb,uuint,vvint,wwint,&
     &x0,y0,z0,xt,yt,zt,nnel,jlev,ttint,ssint,i)
            tsdbt(i,j)=ttint
            ssdbt(i,j)=ssint
          endif
        enddo !j=kbs(i),kfs(i)
478   continue !i=1,ns

      if(nscreen.eq.1) write(*,*)'done 2nd btrack sides'
      write(16,*)'done 2nd btrack sides'

!...  From nodes at half levels
!...
      loop5: do i=1,np

!	Find a wet element to start from
      	index=0
      	do k=1,nne(i)
          ie=ine(i,k)
          if(kfe(ie)/=0) then
            index=1
            ie0=ie
            exit
          endif
        enddo !k

   	if(index==0) then
          do j=kbp(i),kfp(i)
            tndbt(i,j)=tnd(i,j)
            sndbt(i,j)=snd(i,j)
          enddo !j
          cycle loop5
        endif

        do j=kbp(i),kfp(i) !wet node
!         Initialize (x0,y0,z0),nnel and vel.
!	  Caution! nnel must be initialized inside this loop as it is updated inside.
          nnel=ie0
          jlev=j
          x0=x(i)
          y0=y(i)
          if(j.eq.kfe(nnel)) then
            z0=zmsl+eta1(nnel)-dc(nnel,j)/2
          else
            z0=ztot(j)-dc(nnel,j)/2
          endif
!	  Scopes for ww2 etc. have been extended
          uuint=(uu2(i,j)+uu2(i,j-1))/2
          vvint=(vv2(i,j)+vv2(i,j-1))/2
          wwint=(ww2(i,j)+ww2(i,j-1))/2
          vmag=dsqrt(uuint**2+vvint**2+wwint**2)

          if(vmag.le.1.e-4) then
!	  No activity
            tndbt(i,j)=tnd(i,j)
            sndbt(i,j)=snd(i,j)
          else !do backtrack
            ndels=ndelt(i,j)
      	    dtb=dt/ndels !sub-step in backtracking
            call btrack(2,nsplit*ndels,dtb,uuint,vvint,wwint,&
     &x0,y0,z0,xt,yt,zt,nnel,jlev,ttint,ssint,i)
            tndbt(i,j)=ttint
            sndbt(i,j)=ssint
          endif
        enddo !j=kbp(i),kfp(i)
      end do loop5 !i=1,np

      en=0 !secnds(st)
      if(nscreen.eq.1) write(*,*)'2nd btrack took',en,'seconds...'
      write(16,*)'2nd btrack took',en,'seconds...'

!
!*******************************************************************
!								   *
!		     Transport equation				   *
!								   *
!*******************************************************************
!

!... At nodes
!...
!... R.H.S. of transport eq.
!... This part is computed first in order not to interfere with
!... ssd,tsd at new time level.
      do i=1,np
        do j=kbp(i),kfp(i) !wet nodes
          trhs(i,j)=tndbt(i,j)*dzp(i,j)
          srhs(i,j)=sndbt(i,j)*dzp(i,j)

!	  Heat conservation
          if(ihconsv.ne.0) then
            if(j.eq.kfp(i)) then
              zup=zmsl+peta(i)
            else
              zup=ztot(j)
            endif
            zdown=zup-dzp(i,j)
            dp1=dmin1(zmsl+peta(i)-zup,1.d2) !to prevent underflow
            dp2=dmin1(zmsl+peta(i)-zdown,1.d2) !to prevent underflow
            if(dp1.lt.0.or.dp2.lt.0) write(12,*)'Depth<0 in nodal transport',i,dp1,dp2
            sradp1=srad(i)*(0.8*dexp(-dp1/0.9)+0.2*dexp(-dp1/2.1))/rho0/shw
            if(j.eq.kbp(i)) then
              sradp2=0
            else
              sradp2=srad(i)*(0.8*dexp(-dp2/0.9)+0.2*dexp(-dp2/2.1))/rho0/shw
            endif
            trhs(i,j)=trhs(i,j)+nsplit*dt*(sradp1-sradp2)
            if(j.eq.kfp(i)) trhs(i,j)=trhs(i,j)+nsplit*dt*sflux(i)/rho0/shw
          endif !ihconsv.ne.0

        enddo !j=kbp(i),kfp(i)
      enddo !i=1,np


!...
!... solve the transport eqs. and find new temp. and salt
!... 
      do 469 j=1,np
        if(kfp(j).eq.0) go to 469

        nv=kfp(j)-kbp(j)+1

!	Coefficient matrices
!	Middle rows
        do k=kbp(j)+1,kfp(j)-1
          itmp=kfp(j)-k+1 !row #
          alow(itmp)=-tdiffp(j,k)*nsplit*dt/dzphalf(j,k)
          cupp(itmp)=-tdiffp(j,k-1)*nsplit*dt/dzphalf(j,k-1)
          bdia(itmp)=dzp(j,k)-alow(itmp)-cupp(itmp)
        enddo !k

!	1st and last rows
        if(kfp(j).eq.kbp(j)) then
          bdia(1)=dzp(j,kfp(j))
        else 
          alow(nv)=-tdiffp(j,kbp(j))*nsplit*dt/dzphalf(j,kbp(j))
          bdia(nv)=dzp(j,kbp(j))-alow(nv)
          cupp(1)=-tdiffp(j,kfp(j)-1)*nsplit*dt/dzphalf(j,kfp(j)-1)
          bdia(1)=dzp(j,kfp(j))-cupp(1)
        endif

!	RHS
        do k=1,nv
          rrhs(k,1)=trhs(j,kfp(j)+1-k)
          rrhs(k,2)=srhs(j,kfp(j)+1-k)
        enddo

!	Soln for temp. and salt at new level
        call tridag(mnv,nv,2,alow,bdia,cupp,rrhs,soln,gam)
        do k=1,nv
          klev=kfp(j)+1-k
          tnd(j,klev)=soln(k,1)
          snd(j,klev)=soln(k,2)
!	Check against validity range
          if(snd(j,klev).lt.saltmin.or.snd(j,klev).gt.saltmax) then
            write(12,*)'Reset nodal salinity:',it,j,klev,snd(j,klev)
            snd(j,klev)=dmax1(saltmin,dmin1(snd(j,klev),saltmax))
          endif
          if(tnd(j,klev).lt.tempmin.or.tnd(j,klev).gt.tempmax) then
            write(12,*)'Reset nodal temp:',it,j,klev,tnd(j,klev)
            tnd(j,klev)=dmax1(tempmin,dmin1(tnd(j,klev),tempmax))
          endif
        enddo !k=1,nv

469   continue !j=1,np

 
!...  impose temperature, salinity boundary forcings on open bnds
!...  This is only legit if the whole vertical column 
!...  is imposed (otherwise b.c. needs to be imposed before solver).
!...
      do i=1,nope
        do j=1,nond(i)
          nd=iond(i,j)
          do k=1,nvrt
            if(itetype(i).eq.1.or.itetype(i).eq.2) then
              tnd(nd,k)=tth(i)
            else if(itetype(i).eq.3) then
              tnd(nd,k)=tem0(nd,k)
            else if(itetype(i).eq.-1) then
!	    o.b.c.
              isd0=0 !flag
              loop6: do l=1,nne(nd)
                ie=ine(nd,l)
                do ll=1,i34(ie)
                  if(isbs(js(ie,ll)).gt.0) then
                    isd0=js(ie,ll)
                    exit loop6
                  endif
                enddo !ll
              end do loop6 !l
   	      if(isd0.eq.0) then
                write(11,*)'Cannot find an o.b.s.',i,j,nd
               stop
              endif
              if(vn2(isd0,k).lt.0) then
	    	if(kfp(nd)==0) then
		  write(11,*)'Dry open bnd node:',nd
		  stop
		endif
	        tnd(nd,k)=tth(i)
	      endif
            endif
            tnd(nd,k)=dmax1(tempmin,dmin1(tnd(nd,k),tempmax))

            if(isatype(i).eq.1.or.isatype(i).eq.2) then
              snd(nd,k)=sth(i)
            else if(isatype(i).eq.3) then
              snd(nd,k)=sal0(nd,k)
            else if(isatype(i).eq.-1) then
!	  o.b.c.
              isd0=0 !flag
              loop7: do l=1,nne(nd)
                ie=ine(nd,l)
                do ll=1,i34(ie)
                  if(isbs(js(ie,ll)).gt.0) then
                    isd0=js(ie,ll)
                    exit loop7
                  endif
                enddo !ll
              end do loop7 !l
   	      if(isd0.eq.0) then
                write(11,*)'Cannot find an o.b.s.',i,j,nd
                stop
              endif
              if(vn2(isd0,k).lt.0) snd(nd,k)=sth(i)
            endif
            snd(nd,k)=dmax1(saltmin,dmin1(snd(nd,k),saltmax))
          enddo !k=1,nvrt
        enddo !j=1,nond(i)
      enddo !i=1,nope

!...  Extend T,S for wet nodes
      do i=1,np
      	if(kfp(i).eq.0) cycle
  	do j=1,kbp(i)-1
    	  tnd(i,j)=tnd(i,kbp(i))
          snd(i,j)=snd(i,kbp(i))
      	enddo
      	do j=kfp(i)+1,nvrt
          tnd(i,j)=tnd(i,kfp(i))
          snd(i,j)=snd(i,kfp(i))
      	enddo
      enddo !i=1,np

!... At sides
!...
!... R.H.S. of transport eq.
!... This part must be computed first in order not to interfere with
!... ssd,tsd at new time level.
      do i=1,ns
        nd1=isidenode(i,1)
        nd2=isidenode(i,2)
        do j=kbs(i),kfs(i) !wet sides
          trhs(i,j)=tsdbt(i,j)*dz(i,j)
          srhs(i,j)=ssdbt(i,j)*dz(i,j)

!	  Heat conservation
          if(ihconsv.ne.0) then
            if(is(i,2).ne.0) then
              etam=dmax1(eta2(is(i,1)),eta2(is(i,2)))
            else
              etam=eta2(is(i,1))
            endif
            if(j.eq.kfs(i)) then
              zup=zmsl+etam
            else
              zup=ztot(j)
            endif
            zdown=zup-dz(i,j)
            dp1=dmin1(zmsl+etam-zup,1.d2) !to prevent underflow
            dp2=dmin1(zmsl+etam-zdown,1.d2) !to prevent underflow
            if(dp1.lt.0.or.dp2.lt.0) write(12,*)'Depth<0 in side transport',nd1,nd2,dp1,dp2
            srads0=(srad(nd1)+srad(nd2))/2
            srads1=srads0*(0.8*dexp(-dp1/0.9)+0.2*dexp(-dp1/2.1))/rho0/shw
            if(j.eq.kbs(i)) then
              srads2=0
            else
              srads2=srads0*(0.8*dexp(-dp2/0.9)+0.2*dexp(-dp2/2.1))/rho0/shw
            endif
            trhs(i,j)=trhs(i,j)+nsplit*dt*(srads1-srads2)
            if(j.eq.kfs(i)) trhs(i,j)=trhs(i,j)+nsplit*dt*(sflux(nd1)+sflux(nd2))/2/rho0/shw
          endif !ihconsv.ne.0

        enddo !j=kbs(i),kfs(i)
      enddo !i=1,ns

!...
!... solve the transport eqs. and find new temp. and salt
!... 
      do 429 j=1,ns
      	if(kfs(j).eq.0) go to 429

      	nv=kfs(j)-kbs(j)+1

!	Coefficient matrices 
!	Middle rows
        do k=kbs(j)+1,kfs(j)-1
          itmp=kfs(j)-k+1 !row #
          alow(itmp)=-tdiff(j,k)*nsplit*dt/dzhalf(j,k)
          cupp(itmp)=-tdiff(j,k-1)*nsplit*dt/dzhalf(j,k-1)
          bdia(itmp)=dz(j,k)-alow(itmp)-cupp(itmp)
        end do !k

!	1st and last rows
        if(kfs(j).eq.kbs(j)) then
          bdia(1)=dz(j,kfs(j))
        else ! M > m
          alow(nv)=-tdiff(j,kbs(j))*nsplit*dt/dzhalf(j,kbs(j))
          bdia(nv)=dz(j,kbs(j))-alow(nv)
          cupp(1)=-tdiff(j,kfs(j)-1)*nsplit*dt/dzhalf(j,kfs(j)-1)
          bdia(1)=dz(j,kfs(j))-cupp(1)
        endif

!	RHS
      	do k=1,nv
          rrhs(k,1)=trhs(j,kfs(j)+1-k)
          rrhs(k,2)=srhs(j,kfs(j)+1-k)
        enddo

!	Soln for temp. and salt at new level
      	call tridag(mnv,nv,2,alow,bdia,cupp,rrhs,soln,gam)
        do k=1,nv
          klev=kfs(j)+1-k
          tsd(j,klev)=soln(k,1)
          ssd(j,klev)=soln(k,2)
!	Check against validity range
          if(ssd(j,klev).lt.saltmin.or.ssd(j,klev).gt.saltmax) then
!	    write(12,*)'Reset side salinity:',it,j,klev,ssd(j,klev)
            ssd(j,klev)=dmax1(saltmin,dmin1(ssd(j,klev),saltmax))
          endif
          if(tsd(j,klev).lt.tempmin.or.tsd(j,klev).gt.tempmax) then
            write(12,*)'Reset side temp:',it,j,klev,tsd(j,klev)
            tsd(j,klev)=dmax1(tempmin,dmin1(tsd(j,klev),tempmax))
          endif
        enddo !k=1,nv

429   continue !j=1,ns


!...  impose temperature, salinity boundary forcings on open bnds
!...  This is only legit if the whole vertical column 
!...  is imposed (otherwise b.c. needs to be imposed before solver).
!...
      do i=1,nope
        do j=1,nosd(i)
          isd=iosd(i,j)
	  n1=isidenode(isd,1)
	  n2=isidenode(isd,2)
          do k=1,nvrt
            if(itetype(i).eq.1.or.itetype(i).eq.2) then
              tsd(isd,k)=tth(i)
            else if(itetype(i).eq.3) then
              tsd(isd,k)=(tem0(n1,k)+tem0(n2,k))/2
            else if(itetype(i).eq.-1) then
!	o.b.c.
              if(vn2(isd,k).lt.0) tsd(isd,k)=tth(i)
            endif
            tsd(isd,k)=dmax1(tempmin,dmin1(tsd(isd,k),tempmax))

            if(isatype(i).eq.1.or.isatype(i).eq.2) then
              ssd(isd,k)=sth(i)
            else if(isatype(i).eq.3) then
              ssd(isd,k)=(sal0(n1,k)+sal0(n2,k))/2
            else if(isatype(i).eq.-1) then
!	o.b.c.
              if(vn2(isd,k).lt.0) ssd(isd,k)=sth(i)
            endif
            ssd(isd,k)=dmax1(saltmin,dmin1(ssd(isd,k),saltmax))
          enddo !k=1,nvrt
        enddo !j=1,nosd(i)
      enddo !i=1,nope

!...  Extend T,S for wet sides
      do i=1,ns
      	if(kfs(i).eq.0) cycle
  	do j=1,kbs(i)-1
    	  tsd(i,j)=tsd(i,kbs(i))
          ssd(i,j)=ssd(i,kbs(i))
        enddo
        do j=kfs(i)+1,nvrt
          tsd(i,j)=tsd(i,kfs(i))
          ssd(i,j)=ssd(i,kfs(i))
      	enddo
      enddo !i=1,ns

!...  Record vel. at this time step for next internal mode
!...
      do i=1,np
        do k=0,nvrt
          uu1(i,k)=uu2(i,k)
          vv1(i,k)=vv2(i,k)
          ww1(i,k)=ww2(i,k)
        enddo !k
      enddo !i

      if(nscreen.eq.1) write(*,*)'done solving transport equation...'
      write(16,*)'done solving transport equation...'


!... Optional check of heat and salt conservation
      if(ihcheck/=0) then
        open(13,file='hcheck.gr3',status='old')
        read(13,*)
        read(13,*) !ne,np
        do i=1,np
      	  read(13,*)j,xtmp,ytmp,wild
      	  nwild(i)=wild
        enddo !i
        do i=1,ne
	  nwild2(i)=0
	  do j=1,i34(i)
	    if(nwild(nm(i,j))>nwild2(i)) nwild2(i)=nwild(nm(i,j))
	  enddo
        enddo
        close(13)

! 	initialize variables for sub-domain
	temp_freeze = 273.15

! 	heat storage
! 	(Joules)
        heat_storage2 = 0 !new time step
      	salt_storage2 = 0

! 	surface heat fluxes (positive inwards, only those that act at surface)
! 	(Joules)
      	heat_sfc = 0.0 !during previous time step

! 	net solar radiation at surface (positive inwards)
! 	(Joules)
      	solar_sfc = 0.0 !during previous time step

!	Vertical and horizontal advection (positive in)
	heat_vert=0
	heat_adv=0
	salt_vert=0
	salt_adv=0

      	do i=1,ne
          if(nwild2(i)==0) cycle

	  do j=1,i34(i) !nodes
 	    nd=nm(i,j)
	    heat_sfc=heat_sfc+areas(i)*nsplit*dt*sflux(nd)/i34(i)
	    solar_sfc=solar_sfc+areas(i)*nsplit*dt*srad(nd)/i34(i)
	  enddo !j=1,i34 nodes

!	  Advections
	  htmp=0
	  stmp=0
	  icount=0 !for htmp
	  do j=1,i34(i) !sides
 	    isd=js(i,j)
	    if(kfe(i)/=0.and.kfs(isd)/=0) then
	      icount=icount+1
	      htmp=htmp+rho0*shw*we(i,kfe(i))*nsplit*dt*areas(i)*(tsd(isd,kfs(isd))+temp_freeze)
	      stmp=stmp+we(i,kfe(i))*nsplit*dt*areas(i)*ssd(isd,kfs(isd))
	    endif

	    if(ic3(i,j)/=0.and.nwild2(ic3(i,j))==0) then !interface sides
	      do k=kbs(isd),kfs(isd)
	        heat_adv=heat_adv-rho0*shw*vn2(isd,k)*ssign(i,j)*nsplit*dt*distj(isd)*dz(isd,k)* &
     &(tsd(isd,k)+temp_freeze)
	        salt_adv=salt_adv-vn2(isd,k)*ssign(i,j)*nsplit*dt*distj(isd)*dz(isd,k)*ssd(isd,k)
	      enddo !k
	    endif
	  enddo !j=1,i34 sides
	  if(icount/=0) then
	    heat_vert=heat_vert+htmp/icount
	    salt_vert=salt_vert+stmp/icount
	  endif

!	  Heat content
	  do k=kbe(i),kfe(i)
	    hav=0 !average T
	    sav=0
	    icount=0 !for hav
	    do j=1,i34(i) !side
	      isd=js(i,j)
	      if(kfs(isd)/=0) then
	        icount=icount+1
                hav=hav+(tsd(isd,k)+temp_freeze)
                sav=sav+ssd(isd,k)
	      endif
	    enddo !j
	    if(icount/=0) then
	      hav=hav/icount
	      sav=sav/icount
	    endif
	    heat_storage2=heat_storage2+rho0*shw*areas(i)*dc(i,k)*hav
	    salt_storage2=salt_storage2+areas(i)*dc(i,k)*sav
	  enddo !k
        enddo !i=1,ne

        if(it.eq.iths+1) open(90,file='hconsv.dat')
        if(it.eq.iths+1) open(91,file='sconsv.dat')

        write(90,2319)time,heat_storage2,heat_sfc,solar_sfc,heat_vert,heat_adv, &
     &heat_storage2+heat_vert+heat_adv,heat_sfc+solar_sfc
        write(91,2319)time,salt_storage2,salt_vert,salt_adv, &
     &salt_storage2+salt_vert+salt_adv

        if(nscreen.eq.1) write(*,*)'done checking heat & salt budget...'
        write(16,*)'done checking heat & salt budget...'
      endif !ihcheck.ne.0

!----------------------------------------------------------------------
      endif !(ibc.eq.0.or.ibtp.eq.1).and.mod(it,nsplit).eq.0

!...  Compute density at sides, elements and nodes (half levels)
!...
      call eqstate

      if(nscreen.eq.1) write(*,*)'done computing the densities...'
      write(16,*)'done computing the densities...'

!... Optional computation of fluxes and total volume etc.
      if(iflux.ne.0) then
!--------------------------------------------------
!     Compute total mass etc.
      open(10,file='total.dat')
      tvol=0 !total volume
      tmass=0 !total mass
      tpe=0 !total potential energy
      tkne=0 !total kinetic energy
      do i=1,ne
        do j=kbe(i),kfe(i) !wet elements
          tvol=tvol+areas(i)*dc(i,j)  
          tmass=tmass+erho(i,j)*areas(i)*dc(i,j)  
          if(j.eq.kfe(i)) then
            zup=zmsl+eta2(i)
          else
            zup=ztot(j)
          endif
          tpe=tpe+0.5*erho(i,j)*g*areas(i)*dc(i,j)*(2*zup-dc(i,j))

	  uc=0 !average u
	  vc=0 !average u
	  do l=1,i34(i)
	    nd=nm(i,l)
   	    uc=uc+(uu2(nd,j)+uu2(nd,j-1))/2/i34(i)
   	    vc=vc+(vv2(nd,j)+vv2(nd,j-1))/2/i34(i)
	  enddo !l
          wc=(we(i,j)+we(i,j-1))/2 !we has been extended
          tkne=tkne+erho(i,j)*areas(i)*dc(i,j)*(uc**2+vc**2+wc**2)
        enddo !j=kbe(i),kfe(i)
      enddo !i=1,ne
      write(10,*)time/86400,tvol,tmass,tpe,tkne

!	Fluxes
      open(13,file='fluxflag.gr3',status='old')
      read(13,*)
      read(13,*)ntmp,npflag
      if(npflag.ne.np) then
        write(11,*)'# of pts in fluxflag should = np'
       	stop
      endif
      do i=1,np
        read(13,*)j,xtmp,ytmp,wild
        nwild(i)=wild
      enddo
      close(13)
      do i=1,ne
	nwild2(i)=0
	do j=1,i34(i)
	  nd=nm(i,j)
	  if(nwild2(i)<nwild(nd)) nwild2(i)=nwild(nd)
	enddo !j
      enddo !i

      open(9,file='flux.dat')
      open(8,file='salt.dat')
      tvol12=0 !total volume inside rgns 1 and 2
      fluxbnd=0
      fluxchan=0  !total flux across unnatural bnds
      fluxchan1=0  !flux at bnd 1 (positive outward)
      fluxchan2=0  !flux at bnd 2 (positive outward)
      fluxvert=0
      do i=1,ne
        if(nwild2(i)==0) cycle
        if(kfe(i)/=0) then
          fluxvert=fluxvert+we(i,kfe(i))*areas(i)
        endif
        do j=kbe(i),kfe(i) !wet elements
          tvol12=tvol12+areas(i)*dc(i,j)  
        enddo !j=kbe(i),kfe(i)

        do j=1,i34(i)
          isd=js(i,j)
          if(ic3(i,j)==0) then !bnd side
            do k=kbs(isd),kfs(isd) !wet sides
              fluxbnd=fluxbnd+vn2(isd,k)*distj(isd)*dz(isd,k)
            enddo !k
          else if(nwild2(ic3(i,j))==0) then !channel side
            do k=kbs(isd),kfs(isd) !wet sides
              ftmp=ssign(i,j)*vn2(isd,k)*distj(isd)*dz(isd,k)
              fluxchan=fluxchan+ftmp
              if(nwild2(i)==1) then
                fluxchan1=fluxchan1+ftmp
              else !nwild2(i)=2
                fluxchan2=fluxchan2+ftmp
              endif
            enddo !k
          endif !ic3(i,j)=0
        enddo !j=1,i34
      enddo !i=1,ne

      write(9,*)time/86400,fluxchan1,-fluxchan2,fluxbnd,tvol12,fluxchan,&
     &fluxvert,fluxbnd+fluxchan+fluxvert
      if(nscreen.eq.1) write(*,*)'done computing fluxes...'
      write(16,*)'done computing fluxes...'
!---------------------------------------------------------      
      endif !iflux ne 0
!...  end compute flux balance

!... Convert q2 and xl to nodes for outputting
!...
      if(itur.eq.3) then
      	do i=1,np
          do k=kbp(i),nvrt 
	    if(kfp(i).eq.0) then 
	      q2nd(i,k)=-99
              xlnd(i,k)=-99
	    else
              q2nd(i,k)=0
              xlnd(i,k)=0
              do j=1,nne(i)
                ie=ine(i,j)
                id=iself(i,j)
                do l=i34(ie)-2,i34(ie)-1
                  isd1=js(ie,nx(i34(ie),id,l))
                  if(q2(isd1,k).gt.q2nd(i,k)) q2nd(i,k)=q2(isd1,k)
                  if(xl(isd1,k).gt.xlnd(i,k)) xlnd(i,k)=xl(isd1,k)
                enddo !l
              enddo !j
	    endif
          enddo !k
        enddo !i
      endif !itur=3

!...  output global data 
!...
      do j=1,noutput
        if(iof(j).eq.1.and.mod(it,nspool).eq.0) then
	  if(iwrite.eq.0) then
            write(ichan(j),rec=irec(j)+1) real(time)
            write(ichan(j),rec=irec(j)+2) it
            irec(j)=irec(j)+2
	  else !evm
	    a_4 = transfer(source=real(time),mold=a_4)
            write(ichan(j),"(a4)",advance="no") a_4
            a_4 = transfer(source=it,mold=a_4)
            write(ichan(j),"(a4)",advance="no") a_4
	  endif
          do i=1,np
	    if(iwrite.eq.0) then
              write(ichan(j),rec=irec(j)+i) kfp(i)
	    else !evm
	      a_4 = transfer(source=kfp(i),mold=a_4)
              write(ichan(j),"(a4)",advance="no") a_4
	    endif
          enddo !i
          if(iwrite.eq.0) irec(j)=irec(j)+np

          do i=1,np
            if(j.le.10) then
              floatout=0 !in case some are not defined
              if(j.eq.1) then
                floatout=peta(i)
              else if(j.eq.2.and.ihconsv.ne.0) then
                floatout=pr1(i)
              else if(j.eq.3.and.ihconsv.ne.0) then
                floatout=airt1(i)
              else if(j.eq.4.and.ihconsv.ne.0) then
                floatout=shum1(i)
              else if(j.eq.5.and.ihconsv.ne.0) then
                floatout=srad(i)
              else if(j.eq.6.and.ihconsv.ne.0) then
                floatout=fluxsu(i)
              else if(j.eq.7.and.ihconsv.ne.0) then
                floatout=fluxlu(i)
              else if(j.eq.8.and.ihconsv.ne.0) then
                floatout=hradu(i)
              else if(j.eq.9.and.ihconsv.ne.0) then
                floatout=hradd(i)
              else if(j.eq.10.and.ihconsv.ne.0) then
                floatout=sflux(i)
              endif

	      if(iwrite.eq.0) then
                write(ichan(j),rec=irec(j)+1) floatout
                irec(j)=irec(j)+1
	      else !evm
		a_4 = transfer(source=floatout,mold=a_4)
                write(ichan(j),"(a4)",advance="no") a_4
	      endif
            else if(j.eq.11) then
              if(nws.eq.0) then
                windxp(i)=0
                windyp(i)=0
              endif
	      if(iwrite.eq.0) then
                write(ichan(j),rec=irec(j)+1) real(windxp(i))
                write(ichan(j),rec=irec(j)+2) real(windyp(i))
                irec(j)=irec(j)+2
	      else !evm
		a_4 = transfer(source=real(windxp(i)),mold=a_4)
                write(ichan(j),"(a4)",advance="no") a_4
                a_4 = transfer(source=real(windyp(i)),mold=a_4)
                write(ichan(j),"(a4)",advance="no") a_4
	      endif
            else if(j.eq.12) then
              if(nws.ne.2.or.ihconsv.eq.0) then
                floatout=0
                floatout2=0
	      else
                floatout=-tauxz(i)/rho0
                floatout2=-tauyz(i)/rho0
	      endif
	      if(iwrite.eq.0) then
                write(ichan(j),rec=irec(j)+1) real(floatout)
                write(ichan(j),rec=irec(j)+2) real(floatout2)
                irec(j)=irec(j)+2
	      else !evm
		a_4 = transfer(source=real(floatout),mold=a_4)
                write(ichan(j),"(a4)",advance="no") a_4
                a_4 = transfer(source=real(floatout2),mold=a_4)
                write(ichan(j),"(a4)",advance="no") a_4
	      endif
            else if(j.eq.13) then
              do k=kbp(i),nvrt
	        if(iwrite.eq.0) then
                  write(ichan(j),rec=irec(j)+1) real(uu2(i,k))
                  write(ichan(j),rec=irec(j)+2) real(vv2(i,k))
                  irec(j)=irec(j)+2
	    	else !evm
		  a_4 = transfer(source=real(uu2(i,k)),mold=a_4)
                  write(ichan(j),"(a4)",advance="no") a_4
                  a_4 = transfer(source=real(vv2(i,k)),mold=a_4)
                  write(ichan(j),"(a4)",advance="no") a_4
		endif
              enddo !k
            else !j>13
              do k=kbp(i),nvrt
                floatout=0 !for some undefined variables
                if(j.eq.14) then
                  floatout=ww2(i,k)
                else if(kfp(i).eq.0) then
                  floatout=-99
                else
                  if(j.eq.15) then
                    floatout=tnd(i,k)
                  else if(j.eq.16) then
                    floatout=snd(i,k)
                  else if(j.eq.17) then
                    floatout=prho(i,k)
                  else if(j.eq.18) then
                    floatout=tdiffp(i,k)
                  else if(j.eq.19.and.itur.eq.3) then
                    floatout=q2nd(i,k)
                  else if(j.eq.20.and.itur.eq.3) then
                    floatout=xlnd(i,k)
                  endif
                endif

		if(iwrite.eq.0) then
                  write(ichan(j),rec=irec(j)+1) floatout
                  irec(j)=irec(j)+1
		else !evm
		  a_4 = transfer(source=floatout,mold=a_4)
                  write(ichan(j),"(a4)",advance="no") a_4
		endif
              enddo !k
            endif !j
          enddo !i=1,np

          if(nscreen.eq.1) write(*,'(a48)')'done outputting '//variable_nm(j)
          write(16,'(a48)')'done outputting '//variable_nm(j)
        endif !iof(j).eq.1.and.mod(it,nspool).eq.0
      enddo !j=1,noutput

!...  Test output 
      if(noutgm.eq.1.and.mod(it,nspool).eq.0) then
	if(iwrite.eq.0) then
          write(59,rec=igmp+1) real(time)
          write(59,rec=igmp+2) it
          igmp=igmp+2
	else !evm
	  a_4 = transfer(source=real(time),mold=a_4)
          write(59,"(a4)",advance="no") a_4
          a_4 = transfer(source=it,mold=a_4)
          write(59,"(a4)",advance="no") a_4
	endif
        do i=1,ns
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          floatout=0
  	  do j=kbs(i),kfs(i) !wet sides
            uuc=(uu2(n1,j-1)+uu2(n1,j)+uu2(n2,j-1)+uu2(n2,j))/4
            vvc=(vv2(n1,j-1)+vv2(n1,j)+vv2(n2,j-1)+vv2(n2,j))/4
            vnorm=uuc*snx(i)+vvc*sny(i)
            vtan=-uuc*sny(i)+vvc*snx(i)
            diff=dsqrt((vnorm-vn2(i,j))**2+(vtan-vt2(i,j))**2)
            if(diff.gt.floatout) floatout=diff
          enddo !j
	  if(iwrite.eq.0) then
            write(59,rec=igmp+1) floatout
            igmp=igmp+1
	  else
	    a_4 = transfer(source=floatout,mold=a_4)
            write(59,"(a4)",advance="no") a_4
	  endif
        enddo !i
        if(nscreen.eq.1) write(*,*)'done writing test results...'
        write(16,*)'done writing test results...'
      endif !noutgm.eq.1.and.mod(it,nspool).eq.0

!... serg
!...      call da_compute_ugsc !compute uv at side centers in global coordinates
      call da_gout(real(time),it)
!... serg

!...  Open new file
      if(it.ge.ifile*ihfskip) then
	ifile=ifile+1
	write(ifile_char,'(i12)') ifile
	do i=1,noutput
	  close(ichan(i))
	enddo
	call header
!... serg
        do i=1,da_noutput
          close(da_ichan(i))
        enddo
        call da_header
!..  serg
      endif

!...  End output section

!...  write out hot start information if nhstar=1 and at correct time step
!...  note: the hot start files use a record length of 8 on both 32 bit
!...  workstations and the 64 bit cray.
!...
      if(nhstar.eq.1.and.mod(it,ihfskip).eq.0) then
!	Convert it to a string
        write(it_char,'(i12)')it
!... the if condition on treating of da_ihot_type is added by serg
        if (da_ihot_type.eq.1) then !serg
          open(36,file=it_char//'_hotstart',access='direct',recl=ihot_len)
          if(itur.eq.3) then
            write(36,rec=1)time,it,(eta1(i),eta2(i), &
             &(we(i,j),j=0,nvrt),i=1,ne), &
             &((vn2(i,j),vt2(i,j),j=0,nvrt),(tsd(i,j),ssd(i,j),j=1,nvrt),i=1,ns) &
             &,(peta(i),ibad(i),(uu1(i,j),vv1(i,j),ww1(i,j),0,j=0,nvrt), &
             &(tnd(i,j),snd(i,j),tem0(i,j),sal0(i,j),j=1,nvrt),i=1,np), &
             &((q2(i,j),xl(i,j),j=0,nvrt),i=1,ns), &
             &ifile,ifile_char
          else
            write(36,rec=1)time,it,(eta1(i),eta2(i), &
             &(we(i,j),j=0,nvrt),i=1,ne), &
             &((vn2(i,j),vt2(i,j),j=0,nvrt),(tsd(i,j),ssd(i,j),j=1,nvrt),i=1,ns) &
             &,(peta(i),ibad(i),(uu1(i,j),vv1(i,j),ww1(i,j),0,j=0,nvrt), &
             &(tnd(i,j),snd(i,j),tem0(i,j),sal0(i,j),j=1,nvrt),i=1,np), &
             &ifile,ifile_char
          endif
        else if(da_ihot_type.eq.2) then
          call da_hot_out(time,it,itur)
        endif !(da_ihot_type.eq.1) 

        close(36)

        if(nscreen.eq.1) write(*,*) 'Hot start written',it,time
        write(16,*) 'Hot start written',it,time
      endif !end writing hot start file

      if(nscreen.eq.1) write(*,*)'Time step=', it,' Time=',time
      write(16,*)'Time step=', it,' Time=',time


!...  end time stepping loop
!...
      enddo !it=iths,nt


      call date_and_time(date,timestamp)
      if(nscreen.eq.1) write(*,*)'Run completes successfully at ',date,timestamp
      write(16,*)'Run completes successfully at ',date,timestamp

3220  FORMAT(1X,A32,2X,A24,2X,A24)
3645  FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
2120  FORMAT(2X,E20.10,5X,I10)
2453  FORMAT(2(2X,I8),2(2X,E15.8))
2319  FORMAT(10(1x,e23.16))
!                                                                             *
!                                                                             *
!******************************************************************************
!                                                                             *
!                                                                             *
!...  end of program
!...
      stop
      end

!
!********************************************************************************
!	This program solves a tridiagonal system. It was adapted		*
!	from "Numerical recipes in FORTRAN (pp.43 ). 				*
!										*
!	a,b,c,r: input vectors and are not modified by this program.		*
!	       b is the main diagonal, a below and c above. a(1) and c(n) are	*
!	       not used. r is the r.h.s.                        		*
!	nmax: dimension of a etc in the driving routine.			*
!	n: actual rank of the system.						*
!	nc: input; actual # of columns of rhs (<=5). 				*
!	u: output with nc columns (depending on input nc).			*
!	gam: a working array.							*
!********************************************************************************
!

      subroutine tridag(nmax,n,nc,a,b,c,r,u,gam)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)

      integer, intent(in) :: nmax,n,nc
      real(kind=dbl_kind1), dimension(nmax), intent(in) :: a,b,c
      real(kind=dbl_kind1), dimension(nmax,5), intent(in) :: r
      real(kind=dbl_kind1), dimension(nmax), intent(out) :: gam
      real(kind=dbl_kind1), dimension(nmax,5), intent(out) :: u
!      dimension a(nmax),b(nmax),c(nmax),r(nmax,5),u(nmax,5),gam(nmax)

      if(n.lt.1) then
        write(11,*)'Tridag assumes n>=1'
      	stop
      endif
      if(nc.gt.5) then
        write(11,*)'Increase # of columns in tridag'
      	stop
      endif
      if(b(1).eq.0.) then
      	write(11,*)'tridag:  b(1)=0'
      	stop
      endif

      bet=b(1)
      do i=1,nc
        u(1,i)=r(1,i)/bet
      enddo

      do j=2,n
      	gam(j)=c(j-1)/bet
      	bet=b(j)-a(j)*gam(j)
      	if(bet.eq.0) then
          write(11,*)'tridag failed'
          stop
      	endif
        do i=1,nc
          u(j,i)=(r(j,i)-a(j)*u(j-1,i))/bet
        enddo !i
      enddo !j

!     Backsubstitution
      do j=n-1,1,-1
      	do i=1,nc
  	  u(j,i)=u(j,i)-gam(j+1)*u(j+1,i)
      	enddo
      enddo
      
      return
      end

!
!***************************************************************************
!									   *
!     Solve for the density using temperature and salinity along sides.    *
!     From Pond and Pickard's book.
!     validity region: T: [0,40], S: [0:42]				   *
!									   *
!***************************************************************************
!   
      subroutine eqstate
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      logical,save :: first_call=.true.

!...  At nodes
      do 400 i=1,np
        if(kfp(i).eq.0) go to 400 !dry node not updated

        do k=kbp(i),kfp(i)
   	  ttmp=tnd(i,k)
          stmp=snd(i,k)
          if(ttmp.lt.0.or.ttmp.gt.40.or.stmp.lt.0.or.stmp.gt.42) then
            write(11,*)'Invalid temp. or salinity for density'
            write(11,*)ttmp,stmp,i,k,kfp(i)
            stop
          endif

!	Density at one standard atmosphere
          prho(i,k)=1000-0.157406+6.793952E-2*ttmp-9.095290E-3*ttmp**2 &
     &+1.001685E-4*ttmp**3-1.120083E-6*ttmp**4+6.536332E-9*ttmp**5+ &
     &stmp*(0.824493-4.0899E-3*ttmp+&
     &7.6438E-5*ttmp**2-8.2467E-7*ttmp**3+5.3875E-9*ttmp**4)+&
     &dsqrt(stmp)**3*(-5.72466E-3+1.0227E-4*ttmp-1.6546E-6*ttmp**2)+&
     &4.8314E-4*stmp**2

!	Account for pressure for depth > 100m
          if(.not.first_call.and.zmsl-ztot(k).gt.100) then
	    pres=0
	    do kk=k,kfp(i)
	      ptmp=9.81*prho(i,kk)*dzp(i,kk)
	      if(kk.eq.k) then
		pres=pres+1.e-5*ptmp/2
	      else
		pres=pres+1.e-5*ptmp
	      endif
	    enddo !kk 

!	    Secant bulk modulus
	    sbm=19652.21+148.4206*ttmp-2.327105*ttmp**2+1.360477e-2*ttmp**3-5.155288e-5*ttmp**4+ &
     &pres*(3.239908+1.43713e-3*ttmp+1.16092e-4*ttmp**2-5.77905e-7*ttmp**3)+ &
     &pres**2*(8.50935e-5-6.12293e-6*ttmp+5.2787e-8*ttmp**2)+ &
     &stmp*(54.6746-0.603459*ttmp+1.09987e-2*ttmp**2-6.1670e-5*ttmp**3)+ &
     &dsqrt(stmp)**3*(7.944e-2+1.6483e-2*ttmp-5.3009e-4*ttmp**2)+ &
     &pres*stmp*(2.2838e-3-1.0981e-5*ttmp-1.6078e-6*ttmp**2)+1.91075e-4*pres*dsqrt(stmp)**3- &
     &9.9348e-7*pres**2*stmp+2.0816e-8*ttmp*pres**2*stmp+9.1697e-10*ttmp**2*pres**2*stmp
	    ratio=1/(1-pres/sbm)
	    prho(i,k)=prho(i,k)*ratio
	  endif !.not.first_call.and.zmsl-ztot(k).gt.100

!...	  Fake eqstate option
	  if(ieqstate==1) prho(i,k)=1000+0.8*stmp

          if(prho(i,k).lt.980) then
            write(11,*)'Weird density at node',i,k,prho(i,k),ttmp,stmp
            stop
          endif
        enddo !k=kbp(i),kfp(i)
!	Extend to account for level changes
        do k=1,kbp(i)-1
          prho(i,k)=prho(i,kbp(i))
        enddo
        do k=kfp(i)+1,nvrt
          prho(i,k)=prho(i,kfp(i))
        enddo
400   continue !i=1,np

!...  At sides
      do 445 i=1,ns
!	Dry sides will not be updated
      	if(kfs(i).eq.0) go to 445
        do k=kbs(i),kfs(i)
   	  ttmp=tsd(i,k)
          stmp=ssd(i,k)
          if(ttmp.lt.0.or.ttmp.gt.40.or.stmp.lt.0.or.stmp.gt.42) then
            write(11,*)'Invalid side temp. or salinity for density'
            write(11,*)ttmp,stmp,i,k,kfs(i)
            stop
          endif

!	Density at one standard atmosphere
          srho(i,k)=1000-0.157406+6.793952E-2*ttmp-9.095290E-3*ttmp**2&
     &+1.001685E-4*ttmp**3-1.120083E-6*ttmp**4+6.536332E-9*ttmp**5+ &
     &stmp*(0.824493-4.0899E-3*ttmp+&
     &7.6438E-5*ttmp**2-8.2467E-7*ttmp**3+5.3875E-9*ttmp**4)+&
     &dsqrt(stmp)**3*(-5.72466E-3+1.0227E-4*ttmp-1.6546E-6*ttmp**2)+&
     &4.8314E-4*stmp**2

!	Account for pressure for depth > 100m
          if(.not.first_call.and.zmsl-ztot(k).gt.100) then
	    pres=0
	    do kk=k,kfs(i)
	      ptmp=9.81*srho(i,kk)*dz(i,kk)
	      if(kk.eq.k) then
		pres=pres+1.e-5*ptmp/2
	      else
		pres=pres+1.e-5*ptmp
	      endif
	    enddo !kk 

!	    Secant bulk modulus
	    sbm=19652.21+148.4206*ttmp-2.327105*ttmp**2+1.360477e-2*ttmp**3-5.155288e-5*ttmp**4+ &
     &pres*(3.239908+1.43713e-3*ttmp+1.16092e-4*ttmp**2-5.77905e-7*ttmp**3)+ &
     &pres**2*(8.50935e-5-6.12293e-6*ttmp+5.2787e-8*ttmp**2)+ &
     &stmp*(54.6746-0.603459*ttmp+1.09987e-2*ttmp**2-6.1670e-5*ttmp**3)+ &
     &dsqrt(stmp)**3*(7.944e-2+1.6483e-2*ttmp-5.3009e-4*ttmp**2)+ &
     &pres*stmp*(2.2838e-3-1.0981e-5*ttmp-1.6078e-6*ttmp**2)+1.91075e-4*pres*dsqrt(stmp)**3- &
     &9.9348e-7*pres**2*stmp+2.0816e-8*ttmp*pres**2*stmp+9.1697e-10*ttmp**2*pres**2*stmp
	    srho(i,k)=srho(i,k)/(1-pres/sbm)
	  endif !.not.first_call.and.zmsl-ztot(k).gt.100

!...	  Fake eqstate option
	  if(ieqstate==1) srho(i,k)=1000+0.8*stmp

          if(srho(i,k).lt.980) then
            write(11,*)'Weird density at side',i,k,srho(i,k),ttmp,stmp
            stop
          endif
        enddo !k=kbs(i),kfs(i)
!	Extend to account for level changes
        do k=1,kbs(i)-1
          srho(i,k)=srho(i,kbs(i))
        enddo
        do k=kfs(i)+1,nvrt
          srho(i,k)=srho(i,kfs(i))
        enddo
445   continue  !i=1,ns

!...  Compute density at elements
      eslp1: do i=1,ne
      	if(kfe(i).eq.0) cycle eslp1
        do k=kbe(i),kfe(i) !wet elements
          erho(i,k)=0
          icount=0
          do j=1,i34(i)
            nd=nm(i,j)
            if(kfp(nd)/=0) then
              icount=icount+1
              erho(i,k)=erho(i,k)+prho(nd,k) !prho extended
            endif
          enddo !j=1,3

          if(icount==0) then
            write(11,*)'Wet elem. with dry nodes',i
            stop
          else
            erho(i,k)=erho(i,k)/icount
          endif
        enddo !k=kbe(i),kfe(i)
!	Extend to account for level changes
        do k=1,kbe(i)-1
	  erho(i,k)=erho(i,kbe(i))
	enddo !k
        do k=kfe(i)+1,nvrt
	  erho(i,k)=erho(i,kfe(i))
	enddo !k
      end do eslp1 !i=1,ne

      first_call=.false.

      return
      end

!
!********************************************************************
!								    *
!	Routine to update level indices and wetting and drying      *
!								    *
!********************************************************************
!
      subroutine levels(iths,it,itur)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      integer, intent(in) :: iths,it,itur
      integer :: kfs_new(mns),kfp_new(mnp),kfe_new(mne)

!...  Compute new kbe, kfe, dc, dchalf for elements
!...
      do i=1,ne
      	if(it.eq.iths) then !initialize kbe
          kbe(i)=0
          do j=0,nvrt-1
            if(zmsl-dpe(i).ge.ztot(j).and.zmsl-dpe(i).lt.ztot(j+1)) then
              kbe(i)=j+1 !bottom actually at level j, not j+1
              go to 416
            endif
          enddo 
416  	  if(kbe(i).eq.0) then
            write(11,*)'Either zmsl < depth at element',i
            write(11,*)'or land height above highest level'
            write(11,*)'Check fort.19'
            stop
          endif
        endif !it=iths

        htot=dpe(i)+eta2(i)
        if(htot.le.h0) then
          kfe_new(i)=0
        else !not dry
          kfe_new(i)=0
          tmp=zmsl+eta2(i)
          do j=nvrt-1,0,-1
            if(tmp.gt.ztot(j).and.tmp.le.ztot(j+1)) then
              kfe_new(i)=j+1
              dc(i,kfe_new(i))=tmp-ztot(j)
              go to 427
            endif
          enddo !j
427	  if(kfe_new(i).lt.kbe(i)) then
            if(tmp.gt.ztot(nvrt)) then
              write(11,*)'Large elevation at element',i,eta2(i)
            else
              write(11,*)'Impossible 9',i,zmsl+eta2(i),zmsl-dpe(i)
            endif
            stop
          endif

!	Redefine all lower levels to wipe out previous influence
          if(kfe_new(i).eq.kbe(i)) then
            dc(i,kfe_new(i))=htot
            dchalf(i,kfe_new(i))=htot
          else !kfe_new > kbe
!	Next line necessary to ensure bottom level thickness can come back to orginal
            dc(i,kbe(i))=ztot(kbe(i))-(zmsl-dpe(i))
            do j=kbe(i)+1,kfe_new(i)-1
              dc(i,j)=delz(j)
            enddo
            do j=kbe(i),kfe_new(i)-1
              dchalf(i,j)=(dc(i,j+1)+dc(i,j))/2
            enddo
          endif
        endif
             
!	Define w2 (for MYG) for re-wetted elements
        if(it.ne.iths.and.kfe(i).eq.0.and.kfe_new(i).ne.0) then
          do j=0,nvrt
            we(i,j)=0
          enddo !j
        endif !re-wetted
      enddo !i=1,ne

!...  Update f.s. levels
      do i=1,ne
  	kfe(i)=kfe_new(i)
      enddo

!...  Compute new indices for sides
!...
      do i=1,ns
        if(it.eq.iths) then !initialize mmj
          kbs(i)=0
          do j=0,nvrt-1
            if(zmsl-dps(i).ge.ztot(j).and.zmsl-dps(i).lt.ztot(j+1)) then
              kbs(i)=j+1 !bottom actually at level j, not j+1
!              dz(i,kbs(i))=ztot(j+1)-(zmsl-dps(i))
              go to 414
            endif
          enddo 
414  	  if(kbs(i).eq.0) then
            write(11,*)'Either zmsl < depth at side',i
            write(11,*)'or land height above highest level'
            write(11,*)'Check fort.19'
            stop
          endif
        endif !it=iths

      	if(is(i,2).ne.0) then
          etam=dmax1(eta2(is(i,1)),eta2(is(i,2)))
 	else
          etam=eta2(is(i,1))
  	endif
        htot=dps(i)+etam
        if(htot.le.h0) then
          kfs_new(i)=0
        else !not dry
          kfs_new(i)=0
          tmp=zmsl+etam
          do j=nvrt-1,0,-1
            if(tmp.gt.ztot(j).and.tmp.le.ztot(j+1)) then
              kfs_new(i)=j+1
              dz(i,kfs_new(i))=tmp-ztot(j)
              go to 417
            endif
          enddo !j
417	  if(kfs_new(i).lt.kbs(i)) then
            if(tmp.gt.ztot(nvrt)) then
              write(11,*)'Large elevation at side',i,etam
            else
              write(11,*)'Impossible',i,zmsl+etam,zmsl-dps(i)
            endif
            stop
          endif

!	  Redefine all lower levels to wipe out previous influence
          if(kfs_new(i).eq.kbs(i)) then
            dz(i,kfs_new(i))=htot
            dzhalf(i,kfs_new(i)-1)=htot
            dzhalf(i,kfs_new(i))=htot
          else !kfs_new > mmj
!	  Next line necessary to ensure bottom level thickness can come back to orginal
            dz(i,kbs(i))=ztot(kbs(i))-(zmsl-dps(i))
            do j=kbs(i)+1,kfs_new(i)-1
              dz(i,j)=delz(j)
            enddo
            do j=kbs(i),kfs_new(i)-1
              dzhalf(i,j)=(dz(i,j+1)+dz(i,j))/2
            enddo
            dzhalf(i,kfs_new(i))=dz(i,kfs_new(i))
            dzhalf(i,kbs(i)-1)=dz(i,kbs(i))
          endif
        endif
             
!	Define vel., temp. and salt for re-wetted sides
        if(it.ne.iths.and.kfs(i).eq.0.and.kfs_new(i).ne.0) then
          do 451 l=1,nvrt
            vn2(i,l)=0
            vt2(i,l)=0

            tsd(i,l)=0 !initialize
            ssd(i,l)=0
            if(itur.eq.3) then
              q2(i,l)=0
              xl(i,l)=0
            endif
            icount=0
            do j=1,2
              ie=is(i,j)
              if(ie/=0) then
                do k=1,i34(ie)
            	  isd=js(ie,k)
                  if(kfs(isd)/=0) then
                    icount=icount+1
                    lev=min(max(l,kbs(isd)),kfs(isd))
                    tsd(i,l)=tsd(i,l)+tsd(isd,lev)
                    ssd(i,l)=ssd(i,l)+ssd(isd,lev)
       	            if(itur.eq.3) then
                      lev2=min(max(l,kbs(isd)-1),kfs(isd))
                      q2(i,l)=q2(i,l)+q2(isd,lev2)
                      xl(i,l)=xl(i,l)+xl(isd,lev2)
                    endif
                  endif
                enddo !k=1,i34(ie)
              endif !ie.ne.0
            enddo !j=1,2

            if(icount.eq.0) then	
              write(12,*)'Cannot find a wet neighbor side:'
              write(12,*)it,i
	      n1=isidenode(i,1)
	      n2=isidenode(i,2)
              tsd(i,l)=(tem0(n1,l)+tem0(n2,l))/2
              ssd(i,l)=(sal0(n1,l)+sal0(n2,l))/2
              if(itur.eq.3) then
                q2(i,l)=q2min
                xl(i,l)=xlmin1(i)
              endif
            else
              tsd(i,l)=tsd(i,l)/icount
              ssd(i,l)=ssd(i,l)/icount
              if(itur.eq.3) then
                q2(i,l)=q2(i,l)/icount
                xl(i,l)=xl(i,l)/icount
                if(q2(i,l).lt.0.or.xl(i,l).lt.xlmin2(i)) then
                  write(11,*)'Negative q2/xl',it,i,l,q2(i,l),xl(i,l),xlmin2(i)
                  stop
                endif
              endif
            endif
451	  continue !l=1,nvrt

!	  0th level
          vn2(i,0)=0
          vt2(i,0)=0
          if(itur.eq.3) then
            q2(i,0)=q2(i,1)
            xl(i,0)=xl(i,1)
          endif
        endif !re-wetted
      enddo !i=1,ns

      do i=1,ns
  	kfs(i)=kfs_new(i)
      enddo

!...  Compute new nodal levels
!...
      do i=1,np
      	if(it.eq.iths) then !initialize kbp
          kbp(i)=0
          do j=0,nvrt-1
            if(zmsl-dp(i).ge.ztot(j).and.zmsl-dp(i).lt.ztot(j+1)) then
              kbp(i)=j+1 !bottom actually at level j, not j+1
              go to 400
            endif
          enddo 
400  	  if(kbp(i).eq.0) then
            write(11,*)'Either zmsl < depth at side',i
            write(11,*)'or land height above highest level'
            write(11,*)'Check fort.19'
            stop
          endif
        endif !it=iths

        htot=dp(i)+peta(i)
        if(htot.le.h0) then
          kfp_new(i)=0
        else !not dry
          kfp_new(i)=0
          tmp=zmsl+peta(i)
          do j=nvrt-1,0,-1
            if(tmp.gt.ztot(j).and.tmp.le.ztot(j+1)) then
              kfp_new(i)=j+1
              dzp(i,kfp_new(i))=tmp-ztot(j)
              go to 401
            endif
          enddo !j
401	  if(kfp_new(i).lt.kbp(i)) then
            if(tmp.gt.ztot(nvrt)) then
              write(11,*)'Large elevation at node',i,peta(i)
            else
              write(11,*)'Impossible',i,zmsl+peta(i),zmsl-dp(i)
            endif
            stop
          endif

!	  Redefine all lower levels to wipe out previous influence
          if(kfp_new(i).eq.kbp(i)) then
            dzp(i,kfp_new(i))=htot
            dzphalf(i,kfp_new(i))=htot
          else !kfp_new > kbp
!	  Next line necessary to ensure bottom level thickness can come back to orginal
            dzp(i,kbp(i))=ztot(kbp(i))-(zmsl-dp(i))
            do j=kbp(i)+1,kfp_new(i)-1
              dzp(i,j)=delz(j)
            enddo
            do j=kbp(i),kfp_new(i)-1
              dzphalf(i,j)=(dzp(i,j+1)+dzp(i,j))/2
            enddo
!             dzphalf(i,kfp_new(i))=dzp(i,kfp_new(i))
!             dzphalf(i,kbp(i)-1)=dzp(i,kbp(i))
          endif
        endif
             
!	Define temp. and salt for re-wetted nodes
        if(it.ne.iths.and.kfp(i).eq.0.and.kfp_new(i).ne.0) then
          do l=1,nvrt
            tnd(i,l)=0 !initialize
            snd(i,l)=0
            icount=0
            do j=1,nnp(i)
              nd=inp(i,j)
              if(kfp(nd).ne.0) then
                icount=icount+1
                tnd(i,l)=tnd(i,l)+tnd(nd,l) !extended
                snd(i,l)=snd(i,l)+snd(nd,l) !extended
              endif
            enddo !j=1,nnp(i)
            if(icount.eq.0) then	
              write(12,*)'Cannot find a wet neighbor node:'
              write(12,*)it,i
              tnd(i,l)=tem0(i,l)
              snd(i,l)=sal0(i,l)
            else
              tnd(i,l)=tnd(i,l)/icount
              snd(i,l)=snd(i,l)/icount
            endif
          enddo !l=1,nvrt
        endif !re-wetted
      enddo !i=1,np

      do i=1,np
  	kfp(i)=kfp_new(i)
      enddo

!...  Make dry of elements with dry nodes and check element with dry sides
!...  An element will be dry if all nodes or all sides are dry.
      do i=1,ne
	nsum1=0
	nsum2=0
	do j=1,i34(i)
	  nsum1=nsum1+kfp(nm(i,j))
	  nsum2=nsum2+kfs(js(i,j))
	enddo !j
      	if(nsum1==0) kfe(i)=0 
      	if(nsum2==0.and.kfe(i)/=0) then
          write(11,*)'Wet element with 3 dry sides',i
          stop
        endif
      enddo !i

!...  Check dchalf and dc (they are the denominator of terms in vert. momen. eq etc.)
      do i=1,ne
        do j=kbe(i),kfe(i)-1 !wet elements
          if(dchalf(i,j).lt.denom) then
            write(11,*) 'b',i,j,dchalf(i,j),denom
            stop
          endif
        enddo !j=kbe(i),kfe(i)-1

        if(kbe(i).eq.kfe(i).and.dc(i,kfe(i)).lt.h0*(1-1.0e-4)) then !wet elements 
          write(11,*)'Weird single layer el',dc(i,kfe(i)),h0
          stop
        endif
      enddo !i=1,ne

!... Check dzhalf and dz (they are the denominator of terms in momentum eq etc.)
      do i=1,ns
        do j=kbs(i),kfs(i)-1 !wet sides
          if(dzhalf(i,j).lt.denom) then
            write(11,*) i,j,dzhalf(i,j),denom
            stop
          endif
        enddo !j=kbs(i),kfs(i)-1

        if(kbs(i).eq.kfs(i).and.dz(i,kfs(i)).lt.h0*(1-1.0e-4)) then !wet sides
          write(11,*)'Weird single layer',dz(i,kfs(i)),h0
          stop
        endif
      enddo !i=1,ns

!... Check dzphalf and dzp (they are the denominator of terms in transport eq)
      do i=1,np
        do j=kbp(i),kfp(i)-1 !wet nodes
          if(dzphalf(i,j).lt.denom) then
            write(11,*) i,j,dzphalf(i,j),denom
            stop
          endif
        enddo !j=kbp(i),kfp(i)-1

        if(kbp(i).eq.kfp(i).and.dzp(i,kfp(i)).lt.h0*(1-1.0e-4)) then !wet node
          write(11,*)'Weird single layer',dzp(i,kfp(i)),h0
          stop
        endif
      enddo !i=1,np

      return
      end


!
!***************************************************************************
!								     	   *
!     Convert normal vel. to 3D nodal vel. at WHOLE levels.	   	   *
!								     	   *
!***************************************************************************
!
      subroutine nodalvel
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      real(kind=dbl_kind), dimension(0:mnv) :: utmp,vtmp
      real(kind=dbl_kind), dimension(mns,0:mnv) :: ug,vg

!     Side vel. at whole levels
      do i=1,ns
        if(kfs(i).ne.0) then !wet side 
       	  do k=kbs(i),kfs(i) 
            utmp(k)=vn2(i,k)*snx(i)-vt2(i,k)*sny(i)
            vtmp(k)=vn2(i,k)*sny(i)+vt2(i,k)*snx(i)
          enddo !k

          ug(i,kbs(i)-1)=utmp(kbs(i))
          vg(i,kbs(i)-1)=vtmp(kbs(i))
          ug(i,kfs(i))=utmp(kfs(i))
          vg(i,kfs(i))=vtmp(kfs(i))
          do k=kbs(i),kfs(i)-1 !M > m
            ug(i,k)=utmp(k)+dz(i,k)/2/dzhalf(i,k)*(utmp(k+1)-utmp(k))
            vg(i,k)=vtmp(k)+dz(i,k)/2/dzhalf(i,k)*(vtmp(k+1)-vtmp(k))
          enddo !k
        endif !wet side
      enddo !i=1,ns

!     Horizontal vel. at nodes
      do i=1,np
        do k=0,nvrt
          uu2(i,k)=0
          vv2(i,k)=0
          weit=0
          icount=0
          do j=1,nne(i)
            ie=ine(i,j)
            id=iself(i,j)
            do l=i34(ie)-2,i34(ie)-1
              nd=js(ie,nx(i34(ie),id,l))
              if(is(nd,2)==0) then
                fac=2
              else
                fac=1
              endif
              if(kfs(nd)/=0.and.k>=kbs(nd)-1) then
                kin=min(k,kfs(nd))
                uu2(i,k)=uu2(i,k)+ug(nd,kin)/distj(nd)*fac
                vv2(i,k)=vv2(i,k)+vg(nd,kin)/distj(nd)*fac
!	        weit=weit+fac/distj(nd)
                icount=icount+1
              endif
              weit=weit+fac/distj(nd)
            enddo !l
          enddo !j

          if(icount.ne.0) then
            uu2(i,k)=uu2(i,k)/weit
            vv2(i,k)=vv2(i,k)/weit
          endif
        enddo !k=0,nvrt
      enddo !i=1,np

!     w at each node by weighted average
      do i=1,np
        do k=0,nvrt
          ww2(i,k)=0
          twei=0  !sum of weights
          do j=1,nne(i)
            iel=ine(i,j)
            if(kfe(iel).ne.0.and.k.ge.kbe(iel)-1) then
              kin=min(k,kfe(iel))
              twei=twei+areas(iel)
              ww2(i,k)=ww2(i,k)+we(iel,kin)*areas(iel)
            endif
          enddo !j=1,nne(i)

   	  if(twei/=0) ww2(i,k)=ww2(i,k)/twei
        enddo !k=0,nvrt
      enddo !i=1,np

!...  Compute discrepancy between avergaed and elemental vel. vectors 
!      do i=1,np
!	do k=1,nvrt
!	  testa(i,k)=0
!          do j=1,nne(i)
!	    iel=ine(i,j)
!	    index=0
!	    do l=1,3
!	      if(nm(iel,l).eq.i) index=l
!	    enddo !l
!	    if(index.eq.0) then
!	      write(*,*)'Wrong element ball'
!	      stop
!	    endif
!	    testa(i,k)=testa(i,k)+dsqrt((uuf(iel,index,k)-uu2(i,k))**2+
!     +(vvf(iel,index,k)-vv2(i,k))**2)/nne(i)
!	  enddo !j
!	enddo !k
!      enddo !i

      return
      end

!
!***************************************************************************
!									   *
!     			Routine for backtracking. 			   *
!	Input: (x0,y0,z0) and nnel that encloses it, and initial vel.,	   *
!	       initial level (for quicksearch), and a flag indicating      *
!	       1st or 2nd tracking.					   *
!	Output: destination pt (xt,yt,zt), element nnel and level jlev,    *
!	and vel. there (uuint,vvint,wwint), and (S,T) there. 		   *
!									   *
!***************************************************************************
!
      subroutine btrack(imode,ndelt,dtb,uuint,vvint,wwint,x0,y0,z0,&
     &xt,yt,zt,nnel,jlev,ttint,ssint,id0)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      integer, intent(in) :: imode,ndelt,id0
      real(kind=dbl_kind), intent(in) :: dtb
      integer, intent(inout) :: nnel,jlev
      real(kind=dbl_kind), intent(inout) :: uuint,vvint,wwint,x0,y0,z0
      real(kind=dbl_kind), intent(out) :: xt,yt,zt,ttint,ssint
      real(kind=dbl_kind) :: vxl(4,2),vyl(4,2),vzl(4,2),vxn(4),vyn(4),vzn(4)
      real(kind=dbl_kind) :: staint(4),t_xi(4),s_xi(4),sig(4),subrat(4)

      if(imode.ne.1.and.imode.ne.2) then
        write(11,*)'Unknown imode in btrack'
      	stop
      endif

      do idt=1,ndelt
        trat=dble(idt)/ndelt 
        xt=x0-dtb*uuint
        yt=y0-dtb*vvint
        zt=z0-dtb*wwint
      	call quicksearch(1,nnel,jlev,dtb,x0,y0,z0,xt,yt,zt,iflqs1,idt,id0)

!	nnel not dry
!	Interpolate in time
      	do j=1,i34(nnel)
          nd=nm(nnel,j)
          do l=1,2
            lev=jlev+l-2
            if(imode.eq.1) then
              vxl(j,l)=uu2(nd,lev)
              vyl(j,l)=vv2(nd,lev)
              vzl(j,l)=ww2(nd,lev)
            else
              vxl(j,l)=uu2(nd,lev)*(1-trat)+uu1(nd,lev)*trat
              vyl(j,l)=vv2(nd,lev)*(1-trat)+vv1(nd,lev)*trat
              vzl(j,l)=ww2(nd,lev)*(1-trat)+ww1(nd,lev)*trat
            endif
          enddo !l
        enddo !j

!	Interpolate in vertical 
        if(jlev.eq.kfe(nnel)) then
          zup=zmsl+eta1(nnel)
        else
          zup=ztot(jlev)
        endif
        if(zt.gt.zup) then
          write(11,*)'Vertical impossible',zt,zup
          stop
        endif
        zrat=(zup-zt)/dc(nnel,jlev)
      
        do j=1,i34(nnel)
          vxn(j)=vxl(j,2)*(1-zrat)+vxl(j,1)*zrat
          vyn(j)=vyl(j,2)*(1-zrat)+vyl(j,1)*zrat
          vzn(j)=vzl(j,2)*(1-zrat)+vzl(j,1)*zrat
        enddo !j

!	Interpolate in horizontal
        n1=nm(nnel,1)
        n2=nm(nnel,2)
        n3=nm(nnel,3)
	if(i34(nnel)==3) then
          staint(1)=signa(xt,x(n2),x(n3),yt,y(n2),y(n3))/areas(nnel)
          staint(2)=signa(x(n1),xt,x(n3),y(n1),yt,y(n3))/areas(nnel)
          staint(3)=signa(x(n1),x(n2),xt,y(n1),y(n2),yt)/areas(nnel)
          staint(1)=dmax1(0.0d0,dmin1(1.0d0,staint(1)))
          staint(2)=dmax1(0.0d0,dmin1(1.0d0,staint(2)))
          if(staint(1)+staint(2).gt.1) then
            staint(3)=0
            staint(2)=1-staint(1)
          else
            staint(3)=1-staint(1)-staint(2)
          endif
      
!          prat0=dabs(staint1)+dabs(staint2)+dabs(staint3)
	else !quad
	  n4=nm(nnel,4)
          aa1=signa(x(n1),x(n2),xt,y(n1),y(n2),yt)
          aa2=signa(x(n2),x(n3),xt,y(n2),y(n3),yt)
          aa3=signa(x(n3),x(n4),xt,y(n3),y(n4),yt)
          aa4=signa(x(n4),x(n1),xt,y(n4),y(n1),yt)
	  aa=dabs(aa1)+dabs(aa2)+dabs(aa3)+dabs(aa4)
!	  if(aa1<-small1.or.aa2<-small1.or.aa3<-small1.or.aa4<-small1) then
	  if(dabs(aa-areas(nnel))/areas(nnel)>small1) then
	    write(11,*)'Outside quad',aa1,aa2,aa3,aa4,aa,areas(nnel)
	    stop
	  endif
	  call ibilinear(dble(nnel),x(n1),x(n2),x(n3),x(n4),y(n1),y(n2),y(n3),y(n4),xt,yt,csi,etta,staint)
	endif
	uuint=0; vvint=0; wwint=0
        do j=1,i34(nnel)
	  uuint=uuint+vxn(j)*staint(j)
	  vvint=vvint+vyn(j)*staint(j)
	  wwint=wwint+vzn(j)*staint(j)
	enddo !j

        if(iflqs1.eq.1) exit

        x0=xt
        y0=yt
        z0=zt
      enddo !do idt=1,ndelt

!...  Interpolate S,T at the foot
!...  Causion: be very careful about interpolation as large
!...  errors may occur due to small errors!!!
!...

      if(imode.eq.2) then
!----------------------------------------------------------
!     Interpolate S,T along vertical
      zanchor=zup-dc(nnel,jlev)/2
      if(jlev.eq.kfe(nnel).and.zt.gt.zanchor) then
        k1=jlev
        k2=jlev
      	zrat2=0.5 !arbitrary
      else if(zt.gt.zanchor) then !jlev <= M-1
      	k1=jlev
      	k2=jlev+1
     	zrat2=(zt-zanchor)/dchalf(nnel,jlev)
      else !zt <= zanchor
      	k1=max(jlev-1,kbe(nnel))
      	k2=jlev
      	zrat2=1-(zanchor-zt)/dchalf(nnel,k1)
      endif
      zrat2=dmax1(0.0d0,dmin1(1.0d0,zrat2))

      icount=0
      ta=0
      sa=0
      do i=1,i34(nnel)
       	nd=nm(nnel,i)
      	if(kfp(nd).ne.0) then
          icount=icount+1
          t_xi(i)=tnd(nd,k1)*(1-zrat2)+tnd(nd,k2)*zrat2
          s_xi(i)=snd(nd,k1)*(1-zrat2)+snd(nd,k2)*zrat2
          ta=ta+t_xi(i)
          sa=sa+s_xi(i)
        endif
      enddo !i

      if(icount.eq.0) then
        write(11,*)'All dry'
        stop
      else if(icount.ne.i34(nnel)) then !average
        ta=ta/icount
        sa=sa/icount
        ttint=ta
        ssint=sa
      else !Interpolate S,T along horizontal
        ta=ta/icount
        sa=sa/icount
	ttint=0
        ssint=0
	do i=1,i34(nnel)
          ttint=ttint+t_xi(i)*staint(i)
          ssint=ssint+s_xi(i)*staint(i)
	enddo !i
      endif

!      go to 412
!...  Attempt to interpolate inside sub-element
!...
      index=0
      index2=0
      if(i34(nnel)==3) then
        do i=1,4
      	  if(i.le.3) then
            n1=nm(nnel,i)
            n2=js(nnel,nx(i34(nnel),i,2))
            n3=js(nnel,nx(i34(nnel),i,1))
            aa1=signa(xt,xcj(n2),xcj(n3),yt,ycj(n2),ycj(n3))
            aa2=signa(x(n1),xt,xcj(n3),y(n1),yt,ycj(n3))
            aa3=signa(x(n1),xcj(n2),xt,y(n1),ycj(n2),yt)
            aa=dabs(aa1)+dabs(aa2)+dabs(aa3)
            subrat(i)=dabs(aa-areas(nnel)/4)*4/areas(nnel)
            if(subrat(i).lt.100*small1) then
              index=1
              sig(1)=aa1*4/areas(nnel)
              sig(2)=aa2*4/areas(nnel)
              sig(1)=dmax1(0.0d0,dmin1(1.0d0,sig(1)))
              sig(2)=dmax1(0.0d0,dmin1(1.0d0,sig(2)))
              if(sig(1)+sig(2).gt.1) then
                sig(3)=0
                sig(2)=1-sig(1)
              else
                sig(3)=1-sig(1)-sig(2)
              endif

              if(kfp(n1).ne.0.and.kfs(n2).ne.0.and.kfs(n3).ne.0) then
                index2=1
                t_xi(1)=tnd(n1,k1)*(1-zrat2)+tnd(n1,k2)*zrat2
                t_xi(2)=tsd(n2,k1)*(1-zrat2)+tsd(n2,k2)*zrat2
                t_xi(3)=tsd(n3,k1)*(1-zrat2)+tsd(n3,k2)*zrat2
                s_xi(1)=snd(n1,k1)*(1-zrat2)+snd(n1,k2)*zrat2
                s_xi(2)=ssd(n2,k1)*(1-zrat2)+ssd(n2,k2)*zrat2
                s_xi(3)=ssd(n3,k1)*(1-zrat2)+ssd(n3,k2)*zrat2
                ttint=t_xi(1)*sig(1)+t_xi(2)*sig(2)+t_xi(3)*sig(3)
                ssint=s_xi(1)*sig(1)+s_xi(2)*sig(2)+s_xi(3)*sig(3)
                exit
              endif
            endif
          else !i=4
            n1=js(nnel,1)
            n2=js(nnel,2)
            n3=js(nnel,3)
            aa1=signa(xt,xcj(n2),xcj(n3),yt,ycj(n2),ycj(n3))
            aa2=signa(xcj(n1),xt,xcj(n3),ycj(n1),yt,ycj(n3))
            aa3=signa(xcj(n1),xcj(n2),xt,ycj(n1),ycj(n2),yt)
            aa=dabs(aa1)+dabs(aa2)+dabs(aa3)
            subrat(i)=dabs(aa-areas(nnel)/4)*4/areas(nnel)
            if(subrat(i).lt.100*small1) then
              index=1
              sig(1)=aa1*4/areas(nnel)
              sig(2)=aa2*4/areas(nnel)
              sig(1)=dmax1(0.0d0,dmin1(1.0d0,sig(1)))
              sig(2)=dmax1(0.0d0,dmin1(1.0d0,sig(2)))
              if(sig(1)+sig(2).gt.1) then
                sig(3)=0
                sig(2)=1-sig(1)
              else
                sig(3)=1-sig(1)-sig(2)
              endif

              if(kfs(n1).ne.0.and.kfs(n2).ne.0.and.kfs(n3).ne.0) then
                index2=1
                t_xi(1)=tsd(n1,k1)*(1-zrat2)+tsd(n1,k2)*zrat2
                t_xi(2)=tsd(n2,k1)*(1-zrat2)+tsd(n2,k2)*zrat2
                t_xi(3)=tsd(n3,k1)*(1-zrat2)+tsd(n3,k2)*zrat2
                s_xi(1)=ssd(n1,k1)*(1-zrat2)+ssd(n1,k2)*zrat2
                s_xi(2)=ssd(n2,k1)*(1-zrat2)+ssd(n2,k2)*zrat2
                s_xi(3)=ssd(n3,k1)*(1-zrat2)+ssd(n3,k2)*zrat2
                ttint=t_xi(1)*sig(1)+t_xi(2)*sig(2)+t_xi(3)*sig(3)
                ssint=s_xi(1)*sig(1)+s_xi(2)*sig(2)+s_xi(3)*sig(3)
                exit
              endif
            endif
        
          endif
        enddo !i=1,4

      else !quad
        do i=1,4
          n1=nm(nnel,i)
	  n2=js(nnel,nx(i34(nnel),i,3))
	  n4=js(nnel,nx(i34(nnel),i,2))
          aa1=signa(x(n1),xcj(n2),xt,y(n1),ycj(n2),yt)
          aa2=signa(xcj(n2),xctr(nnel),xt,ycj(n2),yctr(nnel),yt)
          aa3=signa(xctr(nnel),xcj(n4),xt,yctr(nnel),ycj(n4),yt)
          aa4=signa(xcj(n4),x(n1),xt,ycj(n4),y(n1),yt)
          bb1=signa(x(n1),xcj(n2),xcj(n4),y(n1),ycj(n2),ycj(n4))
          bb2=signa(xcj(n2),xctr(nnel),xcj(n4),ycj(n2),yctr(nnel),ycj(n4))
	  if(bb1<=0.or.bb2<=0) then
	    write(11,*)'Negative sub-element',nnel,bb1,bb2
	    stop
	  endif
	  bb=bb1+bb2
	  aa=dabs(aa1)+dabs(aa2)+dabs(aa3)+dabs(aa4)
!	  if(aa1>=-small1/10.and.aa2>=-small1/10.and.aa3>=-small1/10.and.aa4>=-small1/10) then
	  subrat(i)=dabs(aa-bb)/bb
	  if(subrat(i) < 4*small1) then
            index=1
	    call ibilinear(nnel+0.1d0*i,x(n1),xcj(n2),xctr(nnel),xcj(n4), &
     &y(n1),ycj(n2),yctr(nnel),ycj(n4),xt,yt,csi,etta,staint)
            if(kfp(n1)/=0.and.kfs(n2)/=0.and.kfs(n4)/=0) then 
              index2=1
              t_xi(1)=tnd(n1,k1)*(1-zrat2)+tnd(n1,k2)*zrat2
              t_xi(2)=tsd(n2,k1)*(1-zrat2)+tsd(n2,k2)*zrat2
              t_xi(3)=ta !icount/=0
              t_xi(4)=tsd(n4,k1)*(1-zrat2)+tsd(n4,k2)*zrat2
              s_xi(1)=snd(n1,k1)*(1-zrat2)+snd(n1,k2)*zrat2
              s_xi(2)=ssd(n2,k1)*(1-zrat2)+ssd(n2,k2)*zrat2
              s_xi(3)=sa
              s_xi(4)=ssd(n4,k1)*(1-zrat2)+ssd(n4,k2)*zrat2
	      ttint=0
	      ssint=0
	      do j=1,4
		ttint=ttint+t_xi(j)*staint(j)
		ssint=ssint+s_xi(j)*staint(j)
	      enddo !j
	    endif !kfp(n1)/=0...
          endif !aa1>=0....
	enddo !i=1,4
	
      endif

      if(index.eq.0) then
      	write(11,*)'Not in any sub-element',nnel,i34(nnel),(subrat(i),i=1,4)
        write(11,109)xt,yt
      	stop
      endif

      if(index2.eq.0.and.icount.eq.i34(nnel)) write(12,*)'All sides are not wet',nnel
412   continue
!----------------------------------------------------------
      endif !(imode.eq.2) 

109   format(2(1x,e26.16))

      return
      end

      function lindex(node,ie)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      lindex=0 !error flag
      do j=1,3
        if(node.eq.nm(ie,j)) lindex=j
      enddo
      if(lindex.eq.0) then
        write(11,*)node,' is not in element ',ie
      	stop
      endif
      
      return
      end

      function sums(x1,x2,x3,x4,y1,y2,y3,y4)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z),integer(i-n)

      sums=dabs((x4-x3)*(y2-y3)+(x2-x3)*(y3-y4))/2+&
     &     dabs((x4-x1)*(y3-y1)-(y4-y1)*(x3-x1))/2+&
     &     dabs((y4-y1)*(x2-x1)-(x4-x1)*(y2-y1))/2
      
      return
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z),integer(i-n)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
      
      return
      end

      subroutine header
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      character(len=48) :: variable_out
!      integer :: ivs,i23d

      do i=1,noutput
        if(iwrite.eq.0) irec(i)=0 !record # for binary
        ichan(i)=59+i !output channel #
        if(i.ge.11.and.i.le.13) then
          ivs=2
        else
          ivs=1
        endif
        if(i.lt.13) then
          i23d=2 !2 or 3D
        else
          i23d=3
        endif

        if(iof(i).eq.1) then
          if(iwrite.eq.0) then
	    open(ichan(i),file=ifile_char//'_'//outfile(i),access='direct',recl=nbyte)
	  else !evm
	    open(ichan(i),file=ifile_char//'_'//outfile(i))
	  endif

	  if(iwrite.eq.0) then
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+m) data_format(nbyte*(m-1)+1:nbyte*m)
            enddo
            irec(i)=irec(i)+48/nbyte
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+m) version(nbyte*(m-1)+1:nbyte*m)
            enddo
            irec(i)=irec(i)+48/nbyte
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+m) start_time(nbyte*(m-1)+1:nbyte*m)
            enddo
            irec(i)=irec(i)+48/nbyte
            variable_out=variable_nm(i)
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+m) variable_out(nbyte*(m-1)+1:nbyte*m)
            enddo
            irec(i)=irec(i)+48/nbyte
            variable_out=variable_dim(i)
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+m) variable_out(nbyte*(m-1)+1:nbyte*m)
            enddo
            irec(i)=irec(i)+48/nbyte

            write(ichan(i),rec=irec(i)+1) nrec
            write(ichan(i),rec=irec(i)+2) real(dt*nspool)
            write(ichan(i),rec=irec(i)+3) nspool
            write(ichan(i),rec=irec(i)+4) ivs
            write(ichan(i),rec=irec(i)+5) i23d
            write(ichan(i),rec=irec(i)+6) real(vpos(i))
            irec(i)=irec(i)+6

	  else !evm
	    write(ichan(i),'(a48)',advance="no") data_format
	    write(ichan(i),'(a48)',advance="no") version
	    write(ichan(i),'(a48)',advance="no") start_time
            write(ichan(i),'(a48)',advance="no") variable_nm(i)	
	    write(ichan(i),'(a48)',advance="no") variable_dim(i)
	    a_4 = transfer(source=nrec,mold=a_4)
            write(ichan(i),"(a4)",advance="no") a_4
            a_4 = transfer(source=real(dt*nspool),mold=a_4)
            write(ichan(i),"(a4)",advance="no") a_4
            a_4 = transfer(source=nspool,mold=a_4)
            write(ichan(i),"(a4)",advance="no") a_4
            a_4 = transfer(source=ivs,mold=a_4)
            write(ichan(i),"(a4)",advance="no") a_4
            a_4 = transfer(source=i23d,mold=a_4)
            write(ichan(i),"(a4)",advance="no") a_4
            a_4 = transfer(source=real(vpos(i)),mold=a_4)
            write(ichan(i),"(a4)",advance="no") a_4
	  endif

!	  Vertical grid
	  if(iwrite.eq.0) then
            write(ichan(i),rec=irec(i)+1) real(zmsl)
            write(ichan(i),rec=irec(i)+2) nvrt
            irec(i)=irec(i)+2
	  else !evm
	    a_4 = transfer(source=real(zmsl),mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=nvrt,mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
	  endif
          do k=1,nvrt
	    if(iwrite.eq.0) then
              write(ichan(i),rec=irec(i)+k) real(ztot(k))
	    else !evm
	      a_4 = transfer(source=real(ztot(k)),mold=a_4)
              write(ichan(i),'(a4)',advance="no") a_4
	    endif
          enddo
          if(iwrite.eq.0) irec(i)=irec(i)+nvrt
          irecm=48/nbyte*5+6+2+nvrt !estimates of total record #

!	  Horizontal grid
	  if(iwrite.eq.0) then
            write(ichan(i),rec=irec(i)+1) np
            write(ichan(i),rec=irec(i)+2) ne
            irec(i)=irec(i)+2
	  else !evm
	    a_4 = transfer(source=np,mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=ne,mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
	  endif
          nsumv=0 !sum up all vertical output levels
          do m=1,np
	    if(iwrite.eq.0) then
              write(ichan(i),rec=irec(i)+1)real(x(m))
              write(ichan(i),rec=irec(i)+2)real(y(m))
              write(ichan(i),rec=irec(i)+3)real(dp(m))
	    else !evm
	      a_4 = transfer(source=real(x(m)),mold=a_4)
              write(ichan(i),'(a4)',advance="no") a_4
              a_4 = transfer(source=real(y(m)),mold=a_4)
              write(ichan(i),'(a4)',advance="no") a_4
              a_4 = transfer(source=real(dp(m)),mold=a_4)
              write(ichan(i),'(a4)',advance="no") a_4
	    endif

!	Bottom indices
  	    ibot=0
            do k=0,nvrt-1
              if(zmsl-dp(m).ge.ztot(k).and.zmsl-dp(m).lt.ztot(k+1)) then
       		ibot=k+1
                exit
     	      endif
            enddo !k
   	    if(ibot.eq.0) then
              write(11,*)'Error in initial bottom',m
              stop
            endif

	    if(iwrite.eq.0) then
              write(ichan(i),rec=irec(i)+4)ibot
              irec(i)=irec(i)+4
	    else !evm
	      a_4 = transfer(source=ibot,mold=a_4)
              write(ichan(i),'(a4)',advance="no") a_4
	    endif
            nsumv=nsumv+nvrt-ibot+1
          enddo !m=1,np
          do m=1,ne
	    if(iwrite.eq.0) then
              write(ichan(i),rec=irec(i)+1)i34(m)
              irec(i)=irec(i)+1
	      do mm=1,i34(m)
                write(ichan(i),rec=irec(i)+mm)nm(m,mm)
	      enddo !mm
              irec(i)=irec(i)+i34(m)
	    else !evm
	      a_4 = transfer(source=i34(m),mold=a_4)
              write(ichan(i),'(a4)',advance="no") a_4
	      do mm=1,i34(m)
                a_4 = transfer(source=nm(m,mm),mold=a_4)
                write(ichan(i),'(a4)',advance="no") a_4
	      enddo !mm
	    endif
          enddo !m

!	  Estimate total # of records
          irecm=irecm+4*np+3*ne
          if(i23d.eq.2) then !2D
            if(ivs.eq.1) then
              irecm=irecm+(2+np+np)*nrec
            else
              irecm=irecm+(2+np+np*2)*nrec
            endif
          else !3D
            if(ivs.eq.1) then
              irecm=irecm+(2+np+nsumv)*nrec
            else
              irecm=irecm+(2+np+nsumv*2)*nrec
            endif
          endif
          if(irecm.gt.mirec) then
            write(11,*)'Output file too large',i,irecm
            stop
          endif
          
          if(iwrite.eq.0) then
	    close(ichan(i)) ! do this to flush the write buffer
            open(ichan(i),file=ifile_char//'_'//outfile(i),access='direct',recl=nbyte)
	  endif
        endif !iof(i)=1
      enddo !i=1,noutput

!     Test output
      if(iwrite.eq.0) igmp=0
      if(noutgm.eq.1) then
	if(iwrite.eq.0) then
          open(59,file=ifile_char//'_test.60',access='direct',recl=nbyte)
          igmp=(4+3+3)*8/nbyte
          write(59,rec=igmp+1) nrec
          write(59,rec=igmp+2) ns
          write(59,rec=igmp+3) real(dt*nspool)
          write(59,rec=igmp+4) nspool
          write(59,rec=igmp+5) 1
          igmp=igmp+5
          close(59)
          open(59,file=ifile_char//'_test.60',access='direct',recl=nbyte)
	else !evm
	  open(59,file=ifile_char//'_test.60')
	  a_4 = transfer(source=nrec,mold=a_4)
          write(59,"(a4)",advance="no") a_4
          a_4 = transfer(source=ns,mold=a_4)
          write(59,"(a4)",advance="no") a_4
          a_4 = transfer(source=real(dt*nspool),mold=a_4)
          write(59,"(a4)",advance="no") a_4
          a_4 = transfer(source=nspool,mold=a_4)
          write(59,"(a4)",advance="no") a_4
          a_4 = transfer(source=1,mold=a_4)
          write(59,"(a4)",advance="no") a_4
	endif
      endif

      return
      end



!******************************************************************************
!                                                                             *
!    Transform from lon,lat (lamda,phi) coordinates into CPP coordinates.     *
!    Lon,Lat must be in radians.                                              *
!                                                                             *
!******************************************************************************

      subroutine cpp(x,y,rlambda,phi,rlambda0,phi0)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)

      r=6378206.4
      x=r*(rlambda-rlambda0)*dcos(phi0)
      y=phi*r

      return
      end

!
!********************************************************************************
!										*
!     Straightline search algorithm. Initially nnel is an element that 	        *
!     encompasses (x0,y0). iloc=0: do not nudge initial pt; iloc=1: nudge.	* 
!     Input: iloc,nnel,x0,y0,z0,xt,yt,zt,jlev, time, and vt2, ww2 for 		*
!	abnormal cases;								*
!     Output the updated end pt (xt,yt,zt) (if so), nnel, jlev, a flag nfl.	*
!     Exit btrack if a bnd or dry element is hit and vel. there is small,	* 
!	or death trap is reached.						*
!										*
!********************************************************************************
!
      subroutine quicksearch(iloc,nnel,jlev,time,x0,y0,z0,&
     &xt,yt,zt,nfl,idt,id0)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      integer, intent(in) :: iloc,idt,id0
      real(kind=dbl_kind), intent(in) :: time,x0,y0,z0
      integer, intent(out) :: nfl
      integer, intent(inout) :: nnel,jlev
      real(kind=dbl_kind), intent(inout) :: xt,yt,zt

      if(iloc.gt.1) then
      	write(11,*)'iloc > 1 not implemented'
      	stop
      endif
      if(kfe(nnel).eq.0.or.zmsl-dpe(nnel).gt.zmsl+eta1(nnel)) then
      	write(11,*)'Starting element is dry'
      	stop
      endif

      nfl=0
      trm=time !time remaining

!     Starting element nel
      nel=nnel
      aa=0
      aa1=0
      do i=1,i34(nel)
        n1=nm(nel,i)
        n2=nm(nel,nx(i34(nel),i,1))
        aa=aa+dabs(signa(x(n1),x(n2),x0,y(n1),y(n2),y0))
        aa1=aa1+dabs(signa(x(n1),x(n2),xt,y(n1),y(n2),yt))
      enddo !i
      ae=dabs(aa-areas(nel))/areas(nel)
      if(ae.gt.small1) then
      	write(11,*)'(x0,y0) not in nnel initially',ae,nnel
      	stop
      endif

      ae=dabs(aa1-areas(nel))/areas(nel)
      if(ae.lt.small1) then
 	nnel=nel
      	go to 400
      endif

!     (xt,yt) not in nel, and thus (x0,y0) and (xt,yt) are distinctive
!     An interior pt close to (x0,y0) to prevent underflow for iloc >=1.
      if(iloc.eq.0) then
        xcg=x0
        ycg=y0
      else if(iloc.eq.1) then
        xcg=(1-1.0d-4)*x0+1.0d-4*xctr(nel)
        ycg=(1-1.0d-4)*y0+1.0d-4*yctr(nel)
!	|x0-xcg|/|x0-xctr|=1.0d-3
      endif

      pathl=dsqrt((xt-xcg)**2+(yt-ycg)**2)
      if(xcg.eq.xt.and.ycg.eq.yt.or.pathl.eq.0) then
      	write(11,*)'Zero path',x0,y0,xt,yt,xcg,ycg
     	stop
      endif

!     Starting edge nel_j
      do j=1,i34(nel)
 	jd1=nm(nel,nx(i34(nel),j,1))
      	jd2=nm(nel,nx(i34(nel),j,2))
        call intersect2(xcg,xt,x(jd1),x(jd2),ycg,yt,y(jd1),y(jd2),iflag,xin,yin,tt1,tt2)
        if(iflag.eq.1) then
          nel_j=j
          go to 399
        endif
      enddo !j=1,3
      write(11,*)'Found no intersecting edges'
      stop

399   continue
 
      zin=z0 !intialize
      it=0
      loop4: do
!----------------------------------------------------------------------------------------
      it=it+1
      if(it.gt.1000) then
        write(12,*)'Death trap reached',idt,id0
      	nfl=1
      	xt=xin
      	yt=yin
      	zt=zin
      	nnel=nel
      	exit loop4
      endif
      md1=nm(nel,nx(i34(nel),nel_j,1))
      md2=nm(nel,nx(i34(nel),nel_j,2))
      
!     Compute z position 
      dist=dsqrt((xin-xt)**2+(yin-yt)**2)
      if(dist/pathl.gt.1+1.0d-4) then
      	write(11,*)'Path overshot'
	stop
      endif
      zin=zt-dist/pathl*(zt-zin)
      trm=trm*dist/pathl !time remaining
      
      pathl=dsqrt((xin-xt)**2+(yin-yt)**2)
      if(pathl.eq.0.or.trm.eq.0) then
        write(11,*)'Target reached'
      	stop
      endif

      lit=0 !flag
!     For horizontal exit and dry elements, compute tangential vel.,
!	update target (xt,yt,zt) and continue.
      if(ic3(nel,nel_j).eq.0.or.kfe(ic3(nel,nel_j)).eq.0) then
        lit=1
        isd=js(nel,nel_j)
      	if(isidenode(isd,1)+isidenode(isd,2).ne.md1+md2) then
          write(11,*)'Wrong side'
          stop
        endif

!       Nudge intersect (xin,yin), and update starting pt
        xin=(1-1.0d-4)*xin+1.0d-4*xctr(nel)
        yin=(1-1.0d-4)*yin+1.0d-4*yctr(nel)
        xcg=xin
        ycg=yin

        xvel=-vt2(isd,jlev)*sny(isd)
        yvel=vt2(isd,jlev)*snx(isd)
        zvel=(ww2(md1,jlev)+ww2(md2,jlev))/2
        xt=xin-xvel*trm
        yt=yin-yvel*trm
        zt=zin-zvel*trm
        hvel=dsqrt(xvel**2+yvel**2)
        if(hvel.lt.1.e-4) then
          nfl=1
          xt=xin
          yt=yin
          zt=zin
          nnel=nel
          exit loop4
        endif
        pathl=hvel*trm
      endif !abnormal cases

!     Search for nel's neighbor with edge nel_j, or in abnormal cases, the same
!     nel
      if(lit.eq.0) nel=ic3(nel,nel_j) !next front element
      aa=0
      do i=1,i34(nel)
        k1=nm(nel,i)
        k2=nm(nel,nx(i34(nel),i,1))
        aa=aa+dabs(signa(x(k1),x(k2),xt,y(k1),y(k2),yt))
      enddo !i
      ae=dabs(aa-areas(nel))/areas(nel)
      if(ae.lt.small1) then
      	nnel=nel
  	exit loop4
      endif

!     Next intersecting edge
      do j=1,i34(nel)
         jd1=nm(nel,nx(i34(nel),j,1))
         jd2=nm(nel,nx(i34(nel),j,2))
         if(jd1.eq.md1.and.jd2.eq.md2.or.jd2.eq.md1.and.jd1.eq.md2) cycle
         call intersect2(xcg,xt,x(jd1),x(jd2),ycg,yt,y(jd1),y(jd2),iflag,xin,yin,tt1,tt2)
         if(iflag.eq.1) then
           nel_j=j !next front edge          
           cycle loop4
         endif
      enddo !j
      write(11,*)'Failed to find next edge',lit,xin,yin,xt,yt,nel,&
     &md1,md2,idt,isidenode(id0,1),isidenode(id0,2)
      stop
!----------------------------------------------------------------------------------------
      end do loop4 

 400  continue
!     No vertical exit from domain
!     Check vertical exit
      if(kfe(nnel).eq.0.or.zmsl-dpe(nnel).gt.zmsl+eta1(nnel)) then
        write(11,*)'Ending element is dry'
        stop
      endif
      zt=dmin1(dmax1(zt,zmsl-dpe(nnel)),zmsl+eta1(nnel))
      if(zt.le.ztot(0)) then
        write(11,*)'Ilegal vertical exit 1!'
        stop
      else if(zt.gt.ztot(nvrt)) then
        write(11,*)'Ilegal vertical exit 2!'
        stop
      else
        do k=0,nvrt-1
          if(zt.gt.ztot(k).and.zt.le.ztot(k+1)) then
            jlev=k+1
            exit
          endif
        enddo
      endif
      
!     Bring jlev into local range (due to small inconsistency in defining kfe and kbe)
      jlev=min0(max0(jlev,kbe(nnel)),kfe(nnel))

      return
      end


!
!********************************************************************************
!										*
!     Program to detect if two segments (1,2) and (3,4) have common pts   	*
!     Assumption: the 4 pts are distinctive.					*
!     The eqs. for the 2 lines are: X=X1+(X2-X1)*tt1 and X=X3+(X4-X3)*tt2.	*
!     Output: iflag: 0: no intersection or colinear; 1: exactly 1 intersection.	*
!     If iflag=1, (xin,yin) is the intersection.				*
!										*
!********************************************************************************
!

      subroutine intersect2(x1,x2,x3,x4,y1,y2,y3,y4,iflag,xin,yin,tt1,tt2)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)
      real(kind=dbl_kind1), parameter :: small2=0.0 !small positive number or 0

      real(kind=dbl_kind1), intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4
      integer, intent(out) :: iflag
      real(kind=dbl_kind1), intent(out) :: xin,yin,tt1,tt2

      tt1=-1000
      tt2=-1000
      iflag=0
      delta=(x2-x1)*(y3-y4)-(y2-y1)*(x3-x4)
      delta1=(x3-x1)*(y3-y4)-(y3-y1)*(x3-x4)
      delta2=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

      if(delta.ne.0.0d0) then
        tt1=delta1/delta
        tt2=delta2/delta
        if(tt1.ge.-small2.and.tt1.le.1+small2.and.tt2.ge.-small2.and.tt2.le.1+small2) then
          iflag=1
          xin=x1+(x2-x1)*tt1
          yin=y1+(y2-y1)*tt1
        endif
      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!												!
!	Inverse bilinear mapping for quadrangles						!
!	Convexity of the quad must have been checked, and (x,y) must not be outside the quad.	!
!												!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ibilinear(elem,x1,x2,x3,x4,y1,y2,y3,y4,x,y,xi,eta,shapef)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)
      real(kind=dbl_kind1), parameter:: small3=1.e-5

      real(kind=dbl_kind1), intent(in) :: elem,x1,x2,x3,x4,y1,y2,y3,y4,x,y !elem for debugging only
      real(kind=dbl_kind1), intent(out) :: xi,eta,shapef(4)

      dimension axi(2),aet(2),bxy(2),root_xi(2),root_et(2)

!...  Consts.
      x0=(x1+x2+x3+x4)/4
      y0=(y1+y2+y3+y4)/4
      axi(1)=x2-x1+x3-x4      
      axi(2)=y2-y1+y3-y4      
      aet(1)=x3+x4-x1-x2
      aet(2)=y3+y4-y1-y2
      bxy(1)=x1-x2+x3-x4
      bxy(2)=y1-y2+y3-y4
!      dxi=bxy(2)*axi(1)-bxy(1)*axi(2)
!      deta=bxy(2)*aet(1)-bxy(1)*aet(2)
      dxi=2*((x3-x4)*(y1-y2)-(y3-y4)*(x1-x2))
      deta=2*((x4-x1)*(y3-y2)-(y4-y1)*(x3-x2))

!...  Inverse mapping
      if(dabs(bxy(1))<small3.and.dabs(bxy(2))<small3.or.dabs(dxi)<small3.and.dabs(deta)<small3) then
        icaseno=1      
!	print*, 'Entering case 1'
	dd=axi(1)*aet(2)-axi(2)*aet(1)
        if(dd==0) then
	  write(11,*)'Case 1 error:',dd
	  stop
	endif
 	xi=4*(aet(2)*(x-x0)-aet(1)*(y-y0))/dd
 	eta=4*(axi(1)*(y-y0)-axi(2)*(x-x0))/dd

      else if(dabs(dxi)<small3.and.dabs(deta)>=small3) then   
        icaseno=2      
!	print*, 'Entering case 2'
 	eta=4*(bxy(2)*(x-x0)-bxy(1)*(y-y0))/deta
	dd=(axi(1)+eta*bxy(1))**2+(axi(2)+eta*bxy(2))**2
	if(dd==0) then
	  write(11,*)'Case 2 error:',dd
          stop
	endif
	xi=((4*(x-x0)-eta*aet(1))*(axi(1)+eta*bxy(1))+(4*(y-y0)-eta*aet(2))*(axi(2)+eta*bxy(2)))/dd

      else if(dabs(dxi)>=small3.and.dabs(deta)<small3) then   
        icaseno=3      
!	print*, 'Entering case 3'
	xi=4*(bxy(2)*(x-x0)-bxy(1)*(y-y0))/dxi
	dd=(aet(1)+xi*bxy(1))**2+(aet(2)+xi*bxy(2))**2
	if(dd==0) then
          write(11,*)'Case 3 error:',dd
          stop
        endif
	eta=((4*(x-x0)-xi*axi(1))*(aet(1)+xi*bxy(1))+(4*(y-y0)-xi*axi(2))*(aet(2)+xi*bxy(2)))/dd

      else !General case
        icaseno=4      
!	print*, 'Entering case 4'
	beta=aet(2)*axi(1)-aet(1)*axi(2)-4*(bxy(2)*(x-x0)-bxy(1)*(y-y0))
	gamma=4*(aet(1)*(y-y0)-aet(2)*(x-x0))
	delta=beta*beta-4*gamma*dxi
	if(delta==0) then
	  xi=-beta/2/dxi
	  eta=(4*(bxy(2)*(x-x0)-bxy(1)*(y-y0))-xi*dxi)/deta
	else if(delta>0) then
!	  print*, 'Entering case 4.2'
	  root_xi(1)=(-beta+dsqrt(delta))/2/dxi
	  root_xi(2)=(-beta-dsqrt(delta))/2/dxi
	  icount=0
	  do i=1,2
	    root_et(i)=(4*(bxy(2)*(x-x0)-bxy(1)*(y-y0))-root_xi(i)*dxi)/deta
	    if(dabs(root_xi(i))<=1.1.and.dabs(root_et(i))<=1.1) then
	      xi=root_xi(i)
	      eta=root_et(i)
	      icount=icount+1
	    endif
	  enddo !i
       	  if(icount==2.and.dabs(root_xi(1)-root_xi(2)).lt.small1) then
!	     Do nothing
!	     xi=root_xi(1)
!            eta=root_et(1)
	  else if(icount/=1) then
	    write(11,*)'Abnormal instances',(root_xi(j),root_et(j),j=1,2),icount,elem
	    write(11,'(10e17.10)')x,y,x1,x2,x3,x4,y1,y2,y3,y4
	    write(11,*)dxi,deta,bxy(1),bxy(2)
	    stop
	  endif

	else
	  write(11,*)'No roots',delta,elem
	  stop
	endif
      endif

      if(dabs(xi)>1.1.or.dabs(eta)>1.1) then
	write(11,*)'Out of bound in ibilinear:',xi,eta,elem,icaseno
	write(11,'(2e17.10)')x,y
        stop
      endif

      xi=dmin1(1.d0,dmax1(xi,-1.d0))
      eta=dmin1(1.d0,dmax1(eta,-1.d0))
      shapef(1)=(1-xi)*(1-eta)/4
      shapef(2)=(1+xi)*(1-eta)/4
      shapef(3)=(1+xi)*(1+eta)/4
      shapef(4)=(1-xi)*(1+eta)/4

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!												!
!	Algebraic Stress Models									!
!												!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine asm(g,i,j,vd,td,qd,qd2)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      integer, intent(in) :: i,j
      real(kind=dbl_kind), intent(in) :: g
      real(kind=dbl_kind), intent(out) :: vd,td,qd,qd2

      if(j==kfs(i)) then
	 bvf=0
	 dundz=(vn2(i,j)-vn2(i,j-1))/dzhalf(i,j-1) !dzhalf defined even for M=m
         dutdz=(vt2(i,j)-vt2(i,j-1))/dzhalf(i,j-1)
      else !j <= M-1 and M > m
	 bvf=g/rho0*(srho(i,j+1)-srho(i,j))/dzhalf(i,j)
	 dundz=(vn2(i,j+1)-vn2(i,j))/dzhalf(i,j)
         dutdz=(vt2(i,j+1)-vt2(i,j))/dzhalf(i,j)
      endif
      shear2=dundz**2+dutdz**2 !same form in 2 frames
      Gh=xl(i,j)**2/2/q2(i,j)*bvf
      Gh=dmin1(dmax1(Gh,-0.28d0),0.0233d0)

      if(stab.eq.'GA') then
        sh=0.49393/(1-34.676*Gh)
        sm=(0.39327-3.0858*Gh)/(1-34.676*Gh)/(1-6.1272*Gh)
	cmiu=dsqrt(2.d0)*sm
	cmiup=dsqrt(2.d0)*sh
	cmiu1=dsqrt(2.d0)*0.2 !for k-eq
	cmiu2=dsqrt(2.d0)*0.2 !for psi-eq.
      else if(stab.eq.'KC') then !Kantha and Clayson
	Ghp=(Gh-(Gh-0.02)**2)/(Gh+0.0233-0.04) !smoothing
	sh=0.4939/(1-30.19*Ghp)
      	sm=(0.392+17.07*sh*Ghp)/(1-6.127*Ghp)
	cmiu=dsqrt(2.d0)*sm
	cmiup=dsqrt(2.d0)*sh
	cmiu1=cmiu/schk
	cmiu2=cmiu/schpsi
      else
	write(11,*)'Unknown ASM:',mid
	stop
      endif

      vd=cmiu*xl(i,j)*dsqrt(q2(i,j))
      td=cmiup*xl(i,j)*dsqrt(q2(i,j))
      qd=cmiu1*xl(i,j)*dsqrt(q2(i,j))
      qd2=cmiu2*xl(i,j)*dsqrt(q2(i,j))

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         sergey's subroutings							!
!										!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!...  coppy of a header subrouting for da outputs
      subroutine da_header
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      character(len=48) :: variable_out
!      integer :: ivs,i23d

      do i=1,da_noutput
        if(iwrite.eq.0) da_irec(i)=0 !record # for binary
        da_ichan(i)=99+i !output channel #
        if (i.eq.1) then
          ivs =1
          i23d=2
        else if (i.eq.2) then
          ivs =1
          i23d=3
        else if (i.eq.3) then
          ivs =1
          i23d=3
        else if (i.eq.4) then
          ivs =2
          i23d=3
        endif

        if(da_iof(i).eq.1) then
          if(iwrite.eq.0) then
	    open(da_ichan(i),file=ifile_char//'_'//da_outfile(i),access='direct',recl=nbyte)
	  else !evm
	    open(da_ichan(i),file=ifile_char//'_'//da_outfile(i))
	  endif

	  if(iwrite.eq.0) then
            do m=1,48/nbyte
              write(da_ichan(i),rec=da_irec(i)+m) data_format(nbyte*(m-1)+1:nbyte*m)
            enddo
            da_irec(i)=da_irec(i)+48/nbyte
            do m=1,48/nbyte
              write(da_ichan(i),rec=da_irec(i)+m) version(nbyte*(m-1)+1:nbyte*m)
            enddo
            da_irec(i)=da_irec(i)+48/nbyte
            do m=1,48/nbyte
              write(da_ichan(i),rec=da_irec(i)+m) start_time(nbyte*(m-1)+1:nbyte*m)
            enddo
            da_irec(i)=da_irec(i)+48/nbyte
            variable_out=da_variable_nm(i)
            do m=1,48/nbyte
              write(da_ichan(i),rec=da_irec(i)+m) variable_out(nbyte*(m-1)+1:nbyte*m)
            enddo
            da_irec(i)=da_irec(i)+48/nbyte
            variable_out=da_variable_dim(i)
            do m=1,48/nbyte
              write(da_ichan(i),rec=da_irec(i)+m) variable_out(nbyte*(m-1)+1:nbyte*m)
            enddo
            da_irec(i)=da_irec(i)+48/nbyte

            write(da_ichan(i),rec=da_irec(i)+1) nrec
            write(da_ichan(i),rec=da_irec(i)+2) real(dt*nspool)
            write(da_ichan(i),rec=da_irec(i)+3) nspool
            write(da_ichan(i),rec=da_irec(i)+4) ivs
            write(da_ichan(i),rec=da_irec(i)+5) i23d
            write(da_ichan(i),rec=da_irec(i)+6) real(da_vpos(i))
            da_irec(i)=da_irec(i)+6

	  else !evm 
	    write(da_ichan(i),'(a48)',advance="no") data_format
	    write(da_ichan(i),'(a48)',advance="no") version
	    write(da_ichan(i),'(a48)',advance="no") start_time
            write(da_ichan(i),'(a48)',advance="no") da_variable_nm(i)	
	    write(da_ichan(i),'(a48)',advance="no") da_variable_dim(i)
	    a_4 = transfer(source=nrec,mold=a_4)
            write(da_ichan(i),"(a4)",advance="no") a_4
            a_4 = transfer(source=real(dt*nspool),mold=a_4)
            write(da_ichan(i),"(a4)",advance="no") a_4
            a_4 = transfer(source=nspool,mold=a_4)
            write(da_ichan(i),"(a4)",advance="no") a_4
            a_4 = transfer(source=ivs,mold=a_4)
            write(da_ichan(i),"(a4)",advance="no") a_4
            a_4 = transfer(source=i23d,mold=a_4)
            write(da_ichan(i),"(a4)",advance="no") a_4
            a_4 = transfer(source=real(da_vpos(i)),mold=a_4)
            write(da_ichan(i),"(a4)",advance="no") a_4
	  endif

!	  Vertical grid
	  if(iwrite.eq.0) then
            write(da_ichan(i),rec=da_irec(i)+1) real(zmsl)
            write(da_ichan(i),rec=da_irec(i)+2) nvrt
            da_irec(i)=da_irec(i)+2
	  else !evm
	    a_4 = transfer(source=real(zmsl),mold=a_4)
            write(da_ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=nvrt,mold=a_4)
            write(da_ichan(i),'(a4)',advance="no") a_4
	  endif
          do k=1,nvrt
	    if(iwrite.eq.0) then
              write(da_ichan(i),rec=da_irec(i)+k) real(ztot(k))
	    else !evm
	      a_4 = transfer(source=real(ztot(k)),mold=a_4)
              write(da_ichan(i),'(a4)',advance="no") a_4
	    endif
          enddo
          if(iwrite.eq.0) da_irec(i)=da_irec(i)+nvrt
          irecm=48/nbyte*5+6+2+nvrt !estimates of total record #
          
          
!	Horizontal grid
        if (i.eq.1) then ! writing element centers grid
	  if(iwrite.eq.0) then
            write(da_ichan(i),rec=da_irec(i)+1) ne
            write(da_ichan(i),rec=da_irec(i)+2) 0
            da_irec(i)=da_irec(i)+2
	  else !evm
	    a_4 = transfer(source=ne,mold=a_4)
            write(da_ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=0,mold=a_4)
            write(da_ichan(i),'(a4)',advance="no") a_4
	  endif
          nsumv=0 !sum up all vertical output levels

          do m=1,ne
	    if(iwrite.eq.0) then
              write(da_ichan(i),rec=da_irec(i)+1)real(xctr(m))
              write(da_ichan(i),rec=da_irec(i)+2)real(yctr(m))
              write(da_ichan(i),rec=da_irec(i)+3)real(dpe(m))
	    else !evm
	      a_4 = transfer(source=real(xctr(m)),mold=a_4)
              write(da_ichan(i),'(a4)',advance="no") a_4
              a_4 = transfer(source=real(yctr(m)),mold=a_4)
              write(da_ichan(i),'(a4)',advance="no") a_4
              a_4 = transfer(source=real(dpe(m)),mold=a_4)
              write(da_ichan(i),'(a4)',advance="no") a_4
	    endif

!	Bottom indices
  	    ibot=0
            do k=0,nvrt-1
              if(zmsl-dpe(m).ge.ztot(k).and.zmsl-dpe(m).lt.ztot(k+1)) then
       		ibot=k+1
                exit
     	      endif
            enddo !k
   	    if(ibot.eq.0) then
              write(11,*)'Error in initial bottom',m
              stop
            endif

	    if(iwrite.eq.0) then
              write(da_ichan(i),rec=da_irec(i)+4)ibot
              da_irec(i)=da_irec(i)+4
	    else !evm
	      a_4 = transfer(source=ibot,mold=a_4)
              write(da_ichan(i),'(a4)',advance="no") a_4
	    endif
            nsumv=nsumv+nvrt-ibot+1
          enddo !m=1,ne
        else
!... writnig cide centers grid
         if(iwrite.eq.0) then
            write(da_ichan(i),rec=da_irec(i)+1) ns
            write(da_ichan(i),rec=da_irec(i)+2) 0
            da_irec(i)=da_irec(i)+2
          else !evm
            a_4 = transfer(source=ns,mold=a_4)
            write(da_ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=0,mold=a_4)
            write(da_ichan(i),'(a4)',advance="no") a_4
          endif
          nsumv=0 !sum up all vertical output levels

          do m=1,ns
            if(iwrite.eq.0) then
              write(da_ichan(i),rec=da_irec(i)+1)real(xcj(m))
              write(da_ichan(i),rec=da_irec(i)+2)real(ycj(m))
              write(da_ichan(i),rec=da_irec(i)+3)real(dps(m))
            else !evm
              a_4 = transfer(source=real(xcj(m)),mold=a_4)
              write(da_ichan(i),'(a4)',advance="no") a_4
              a_4 = transfer(source=real(ycj(m)),mold=a_4)
              write(da_ichan(i),'(a4)',advance="no") a_4
              a_4 = transfer(source=real(dps(m)),mold=a_4)
              write(da_ichan(i),'(a4)',advance="no") a_4
            endif

!       Bottom indices
            ibot=0
            do k=0,nvrt-1
              if(zmsl-dps(m).ge.ztot(k).and.zmsl-dps(m).lt.ztot(k+1)) then
                ibot=k+1
                exit
              endif
            enddo !k
            if(ibot.eq.0) then
              write(11,*)'Error in initial bottom',m
              stop
            endif

            if(iwrite.eq.0) then
              write(da_ichan(i),rec=da_irec(i)+4)ibot
              da_irec(i)=da_irec(i)+4
            else !evm
              a_4 = transfer(source=ibot,mold=a_4)
              write(da_ichan(i),'(a4)',advance="no") a_4
            endif
            nsumv=nsumv+nvrt-ibot+1
          enddo !m=1,ns
        endif ! i.eq.1

          if(iwrite.eq.0) then
	    close(da_ichan(i)) ! do this to flush the write buffer
            open(da_ichan(i),file=ifile_char//'_'//da_outfile(i),access='direct',recl=nbyte)
	  endif
        endif !iof(i)=1
      enddo !i=1,noutput

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!...   output da varaibles				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine da_gout(time,it)
    use global
      
!    print *, it, time
      do j=1,da_noutput
        if(da_iof(j).eq.1.and.mod(it,nspool).eq.0) then
	  if(iwrite.eq.0) then
            write(da_ichan(j),rec=da_irec(j)+1) real(time)
            write(da_ichan(j),rec=da_irec(j)+2) it
            da_irec(j)=da_irec(j)+2
	  else !evm
	    a_4 = transfer(source=real(time),mold=a_4)
            write(da_ichan(j),"(a4)",advance="no") a_4
            a_4 = transfer(source=it,mold=a_4)
            write(da_ichan(j),"(a4)",advance="no") a_4
	  endif

         if (j.eq.1) then
           do i=1,ne
	     if(iwrite.eq.0) then
               write(da_ichan(j),rec=da_irec(j)+i) kfe(i)
	     else !evm
	       a_4 = transfer(source=kfe(i),mold=a_4)
               write(da_ichan(j),"(a4)",advance="no") a_4
	     endif
           enddo !i
           if(iwrite.eq.0) da_irec(j)=da_irec(j)+ne

           do i=1,ne
               if(j.eq.1) then
                 floatout=eta2(i)
               endif

	       if(iwrite.eq.0) then
                 write(da_ichan(j),rec=da_irec(j)+1) floatout
                 da_irec(j)=da_irec(j)+1
   	       else !evm
	         a_4 = transfer(source=floatout,mold=a_4)
                 write(da_ichan(j),"(a4)",advance="no") a_4
	       endif

           enddo !i=1,ne

         else if(j.eq.4) then 

           call da_compute_ugsc !compute uv at side centers in global coordinates
           do i=1,ns
             if(iwrite.eq.0) then
               write(da_ichan(j),rec=da_irec(j)+i) kfs(i)
             else !evm
               a_4 = transfer(source=kfs(i),mold=a_4)
               write(da_ichan(j),"(a4)",advance="no") a_4
             endif
           enddo !i
           if(iwrite.eq.0) da_irec(j)=da_irec(j)+ns
           
           do i=1,ns
               do k=kbs(i),nvrt
                 if(iwrite.eq.0) then
                   write(da_ichan(j),rec=da_irec(j)+1) real(da_ugsc(i,k))
                   write(da_ichan(j),rec=da_irec(j)+2) real(da_vgsc(i,k))
                   da_irec(j)=da_irec(j)+2
                 else !evm
                   a_4 = transfer(source=real(vn2(i,k)),mold=a_4)
                   write(da_ichan(j),"(a4)",advance="no") a_4
                   a_4 = transfer(source=real(vt2(i,k)),mold=a_4)
                   write(da_ichan(j),"(a4)",advance="no") a_4
                 endif
               enddo !k
           enddo !i

         else !j~=1.or.4
           do i=1,ns
             if(iwrite.eq.0) then
               write(da_ichan(j),rec=da_irec(j)+i) kfs(i)
             else !evm
               a_4 = transfer(source=kfs(i),mold=a_4)
               write(da_ichan(j),"(a4)",advance="no") a_4
             endif
           enddo !i
           if(iwrite.eq.0) da_irec(j)=da_irec(j)+ns
           
           do i=1,ns
               do k=kbs(i),nvrt
 	         if(j.eq.2) then
            	   floatout=tsd(i,k)
          	 else if(j.eq.3) then
          	   floatout=ssd(i,k)
          	 endif !j.eq.

                 if(iwrite.eq.0) then
                   write(da_ichan(j),rec=da_irec(j)+1) real(floatout)
                   da_irec(j)=da_irec(j)+1
                 else !evm
                   a_4 = transfer(source=real(floatout),mold=a_4)
                   write(da_ichan(j),"(a4)",advance="no") a_4
                 endif
               enddo !k
           enddo !i
         endif !j==1

          if(nscreen.eq.1) write(*,'(a48)')'done outputting '//da_variable_nm(j)
          write(16,'(a48)')'done outputting '//da_variable_nm(j)
        endif !iof(j).eq.1.and.mod(it,nspool).eq.0
      enddo !j=1,noutput
return 
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	convert vn2 vt2 to global coord		!
!	changes glob vars:  da_ugsc da_vgsc	!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine da_compute_ugsc
      use global
        do i=1,ns
          do k=1,nvrt
            da_ugsc(i,k)=vn2(i,k)*snx(i)-vt2(i,k)*sny(i)
            da_vgsc(i,k)=vn2(i,k)*sny(i)+vt2(i,k)*snx(i)
          enddo !k
        enddo !i
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	outputs hotstart in da format
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine da_hot_out(time,it,itur)
      
      use global
      real(kind=dbl_kind), intent(in) :: time
      integer, intent(in) :: it, itur
      character(len=12) :: it_char
!      print *,time,it,itur
       
         call da_compute_ugsc

         ihot_len_da=nbyte*(3+4*ne+2*ne*(nvrt+1)+4*ns*(nvrt+1)+4*ns*nvrt+ &
         &         8*np*nvrt+1)+12
         if(itur.eq.3) ihot_len_da=ihot_len_da+nbyte*4*ns*(nvrt+1)
         write(it_char,'(i12)')it

         open(36,file=it_char//'_hotstart',access='direct',recl=ihot_len_da)
         if(itur.eq.3) then
           write(36,rec=1)time,it,(eta1(i),eta2(i), &
            &(we(i,j),j=0,nvrt),i=1,ne), &
            &((da_ugsc(i,j),da_vgsc(i,j),j=0,nvrt),(tsd(i,j),ssd(i,j),j=1,nvrt),i=1,ns) &
            &,((tnd(i,j),snd(i,j),tem0(i,j),sal0(i,j),j=1,nvrt),i=1,np), &
            &((q2(i,j),xl(i,j),j=0,nvrt),i=1,ns), &
            &ifile,ifile_char
         else
           write(36,rec=1)time,it,(eta1(i),eta2(i), &
            &(we(i,j),j=0,nvrt),i=1,ne), &
            &((da_ugsc(i,j),da_vgsc(i,j),j=0,nvrt),(tsd(i,j),ssd(i,j),j=1,nvrt),i=1,ns) &
            &,((tnd(i,j),snd(i,j),tem0(i,j),sal0(i,j),j=1,nvrt),i=1,np), &
            &ifile,ifile_char
         endif
      end subroutine

